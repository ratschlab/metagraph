#pragma once

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/vector_map.hpp"
#include "common/vector.hpp"
#include "common/logger.hpp"
#include "common/utils/template_utils.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace annot {
namespace binmat {

const std::string kRowDiffAnchorExt = ".anchors";
const std::string kRowDiffForkSuccExt = ".rd_succ";

const size_t RD_PATH_RESERVE_SIZE = 2;


class IRowDiff {
  public:
    virtual ~IRowDiff() {}

    const graph::DBGSuccinct* graph() const { return graph_; }
    void set_graph(const graph::DBGSuccinct *graph) { graph_ = graph; }

  protected:
    const graph::DBGSuccinct *graph_ = nullptr;
};

/**
 * Sparsified representation of the underlying #BinaryMatrix that stores diffs between
 * successive nodes, rather than the full annotation.
 * The successor of a node (that is the node to diff against) is determined by a path in
 * an external graph structure, a #graph::DBGSuccinct, which has the property that rows
 * on a path are likely identical or very similar.
 *
 * RowDiff sparsification can be applied to any BinaryMatrix instance.
 * The row-diff binary matrix is defined by three data structures:
 *   1. #diffs_ the underlying sparsified (diffed) #BinaryMatrix
 *   2. #anchor_ rows marked as anchors are stored in full
 *   3. #graph_ the graph that was used to determine adjacent rows for sparsification
 * Retrieving data from RowDiff requires the associated #graph_. In order to get the
 * annotation for  i-th row, we start traversing the node corresponding to i in #graph_
 * and accumulate the values in #diffs until we hit an anchor node, which is stored in
 * full.
 */
//NOTE: Clang aggressively abuses the clause in the C++ standard (14.7.1/11) that allows
// virtual methods in template classes to not be instantiated if unused and mistakenly
// does not instantiate the virtual methods in this class, so I had to move definitions
// to the header (gcc works fine)
template <class BaseMatrix>
class RowDiff : public IRowDiff, public BinaryMatrix {
  public:
    using anchor_bv_type = bit_vector_small;
    using fork_succ_bv_type = bit_vector_small;

    RowDiff() {}

    RowDiff(const graph::DBGSuccinct *graph, BaseMatrix &&diff)
          : diffs_(std::move(diff)) { graph_ = graph; }

    /**
     * Returns the number of set bits in the row-diff transformed matrix.
     */
    uint64_t num_relations() const override { return diffs_.num_relations(); }
    uint64_t num_columns() const override { return diffs_.num_columns(); }
    uint64_t num_rows() const override { return diffs_.num_rows(); }

    bool get(Row row, Column column) const override;

    /**
     * Returns the given column.
     */
    std::vector<Row> get_column(Column column) const override;

    SetBitPositions get_row(Row row) const override;

    std::vector<SetBitPositions> get_rows(const std::vector<Row> &row_ids) const override;

    bool load(std::istream &f) override;
    void serialize(std::ostream &f) const override;

    void load_fork_succ(const std::string &filename);
    void load_anchor(const std::string &filename);
    const anchor_bv_type& anchor() const { return anchor_; }

    const BaseMatrix& diffs() const { return diffs_; }
    BaseMatrix& diffs() { return diffs_; }

  private:
    static void add_diff(const Vector<uint64_t> &diff, Vector<uint64_t> *row);

    BaseMatrix diffs_;
    anchor_bv_type anchor_;
    fork_succ_bv_type fork_succ_;
};

template <class BaseMatrix>
bool RowDiff<BaseMatrix>::get(Row row, Column column) const {
    assert(graph_ && "graph must be loaded");
    assert(anchor_.size() == diffs_.num_rows() && "anchors must be loaded");
    assert(!fork_succ_.size() || fork_succ_.size() == graph_->num_nodes() + 1);

    SetBitPositions set_bits = get_row(row);
    SetBitPositions::iterator v = std::lower_bound(set_bits.begin(), set_bits.end(), column);
    return v != set_bits.end() && *v == column;
}


/**
 * Returns the given column.
 */
template <class BaseMatrix>
std::vector<BinaryMatrix::Row> RowDiff<BaseMatrix>::get_column(Column column) const {
    assert(graph_ && "graph must be loaded");
    assert(anchor_.size() == diffs_.num_rows() && "anchors must be loaded");
    assert(!fork_succ_.size() || fork_succ_.size() == graph_->num_nodes() + 1);

    // TODO: implement a more efficient algorithm
    std::vector<Row> result;
    for (Row row = 0; row < num_rows(); ++row) {
        if (get(row, column))
            result.push_back(row);
    }
    return result;
}

template <class BaseMatrix>
BinaryMatrix::SetBitPositions RowDiff<BaseMatrix>::get_row(Row row) const {
    assert(graph_ && "graph must be loaded");
    assert(anchor_.size() == diffs_.num_rows() && "anchors must be loaded");
    assert(!fork_succ_.size() || fork_succ_.size() == graph_->num_nodes() + 1);

    Vector<uint64_t> result = diffs_.get_row(row);
    std::sort(result.begin(), result.end());

    auto node = graph::AnnotatedSequenceGraph::anno_to_graph_index(row);

    while (!anchor_[row]) {
        node = graph_->row_diff_successor(node, fork_succ_);
        row = graph::AnnotatedSequenceGraph::graph_to_anno_index(node);

        auto diff_row = diffs_.get_row(row);
        std::sort(diff_row.begin(), diff_row.end());
        add_diff(diff_row, &result);
    }

    return result;
}

template <class BaseMatrix>
std::vector<BinaryMatrix::SetBitPositions>
RowDiff<BaseMatrix>::get_rows(const std::vector<Row> &row_ids) const {
    assert(graph_ && "graph must be loaded");
    assert(anchor_.size() == diffs_.num_rows() && "anchors must be loaded");
    assert(!fork_succ_.size() || fork_succ_.size() == graph_->num_nodes() + 1);

    // diff rows annotating nodes along the row-diff paths
    std::vector<Row> rd_ids;
    rd_ids.reserve(row_ids.size() * RD_PATH_RESERVE_SIZE);

    // map row index to its index in |rd_rows|
    VectorMap<Row, size_t> node_to_rd;
    node_to_rd.reserve(row_ids.size() * RD_PATH_RESERVE_SIZE);

    // Truncated row-diff paths, indexes to |rd_rows|.
    // The last index in each path points to an anchor or to a row which had
    // been reached before, and thus, will be reconstructed before this one.
    std::vector<std::vector<size_t>> rd_paths_trunc(row_ids.size());

    for (size_t i = 0; i < row_ids.size(); ++i) {
        Row row = row_ids[i];

        while (true) {
            auto [it, is_new] = node_to_rd.try_emplace(row, rd_ids.size());
            rd_paths_trunc[i].push_back(it.value());

            // If a node had been reached before, we interrupt the diff path.
            // The annotation for that node will have been reconstructed earlier
            // than for other nodes in this path as well. Thus, we will start
            // reconstruction from that node and don't need its successors.
            if (!is_new)
                break;

            rd_ids.push_back(row);

            if (anchor_[row])
                break;

            auto node = graph::AnnotatedSequenceGraph::anno_to_graph_index(row);
            node = graph_->row_diff_successor(node, fork_succ_);
            row = graph::AnnotatedSequenceGraph::graph_to_anno_index(node);
        }
    }

    node_to_rd = VectorMap<Row, size_t>();

    std::vector<SetBitPositions> rd_rows = diffs_.get_rows(rd_ids);

    rd_ids = std::vector<Row>();

    // reconstruct annotation rows from row-diff
    std::vector<SetBitPositions> rows(row_ids.size());

    for (size_t i = 0; i < row_ids.size(); ++i) {
        SetBitPositions &result = rows[i];
        // propagate back and reconstruct full annotations for predecessors
        for (auto it = rd_paths_trunc[i].rbegin(); it != rd_paths_trunc[i].rend(); ++it) {
            std::sort(rd_rows[*it].begin(), rd_rows[*it].end());
            add_diff(rd_rows[*it], &result);
            // replace diff row with full reconstructed annotation
            rd_rows[*it] = result;
        }
    }

    return rows;
}

template <class BaseMatrix>
bool RowDiff<BaseMatrix>::load(std::istream &f) {
    auto pos = f.tellg();
    std::string version(4, '\0');
    if (f.read(version.data(), 4) && version == "v2.0") {
        if constexpr(!std::is_same_v<BaseMatrix, ColumnMajor>) {
            if (!anchor_.load(f) || !fork_succ_.load(f))
                return false;
        }
    } else {
        // backward compatibility
        f.seekg(pos);
        if constexpr(!std::is_same_v<BaseMatrix, ColumnMajor>) {
            if (!anchor_.load(f))
                return false;

            common::logger->warn(
                "Loading old version of RowDiff without a fork routing bitmap."
                " The last outgoing edges will be used as successors.");
            fork_succ_ = fork_succ_bv_type();
        }
    }
    return diffs_.load(f);
}

template <class BaseMatrix>
void RowDiff<BaseMatrix>::serialize(std::ostream &f) const {
    f.write("v2.0", 4);
    if constexpr(!std::is_same_v<BaseMatrix, ColumnMajor>) {
        anchor_.serialize(f);
        fork_succ_.serialize(f);
    }
    diffs_.serialize(f);
}

template <class BaseMatrix>
void RowDiff<BaseMatrix>::add_diff(const Vector<uint64_t> &diff, Vector<uint64_t> *row) {
    assert(std::is_sorted(row->begin(), row->end()));
    assert(std::is_sorted(diff.begin(), diff.end()));

    if (diff.empty())
        return;

    Vector<uint64_t> result;
    result.reserve(row->size() + diff.size());
    std::set_symmetric_difference(row->begin(), row->end(),
                                  diff.begin(), diff.end(),
                                  std::back_inserter(result));
    row->swap(result);
}

} // namespace binmat
} // namespace annot
} // namespace mtg
