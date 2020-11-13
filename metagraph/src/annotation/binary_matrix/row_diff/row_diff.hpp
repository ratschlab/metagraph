#pragma once

#include <fstream>
#include <string>
#include <vector>

#include <sdsl/enc_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>

#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/logger.hpp"
#include "common/utils/template_utils.hpp"
#include "common/vector.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"

namespace mtg {
namespace annot {
namespace binmat {

const std::string kRowDiffAnchorExt = ".anchors";

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
 *   2. #terminal_ rows marked as terminal are stored in full
 *   3. #graph_ the graph that was used to determine adjacent rows for sparsification
 * Retrieving data from RowDiff requires the associated #graph_. In order to get the
 * annotation for  i-th row, we start traversing the node corresponding to i in #graph_
 * and accumulate the values in #diffs until we hit a terminal node, which is stored in
 * full.
 */
//NOTE: Clang aggressively abuses the clause in the C++ standard (14.7.1/11) that allows
// virtual methods in template classes to not be instantiated if unused and mistakenly
// does not instantiate the virtual methods in this class, so I had to move definitions
// to the header (gcc works fine)
template <class BaseMatrix>
class RowDiff : public BinaryMatrix {
  public:
    using anchor_bv_type = bit_vector_small;

    RowDiff() {}

    RowDiff(const graph::DBGSuccinct *graph, BaseMatrix &&diff)
        : graph_(graph), diffs_(std::move(diff)) {}

    uint64_t num_columns() const override { return diffs_.num_columns(); }

    const graph::DBGSuccinct* graph() const { return graph_; }

    /**
     * Returns the number of set bits in the matrix.
     */
    uint64_t num_relations() const override { return diffs_.num_relations(); }

    uint64_t num_rows() const override { return diffs_.num_rows(); }
    void set_graph(const graph::DBGSuccinct *graph) { graph_ = graph; }

    bool get(Row row, Column column) const override;

    /**
     * Returns the given column.
     */
    std::vector<Row> get_column(Column column) const override;

    SetBitPositions get_row(Row row) const override;

    inline bool load(std::istream &f) override;
    inline void serialize(std::ostream &f) const override;

    void serialize(const std::string &filename) const;
    bool load(const std::string &filename);

    void load_anchor(const std::string& filename);
    const anchor_bv_type& anchor() const { return anchor_; }

    const BaseMatrix& diffs() const { return diffs_; }
    BaseMatrix& diffs() { return diffs_; }

    Vector<uint64_t> get_diff(uint64_t node_id) const { return diffs_.get_row(node_id); }

  private:
    static void merge(Vector<uint64_t> *result, const Vector<uint64_t> &diff2);

    const graph::DBGSuccinct *graph_ = nullptr;

    BaseMatrix diffs_;
    anchor_bv_type anchor_;
};

template <class BaseMatrix>
bool RowDiff<BaseMatrix>::get(Row row, Column column) const {
    assert(anchor_.size() == diffs_.num_rows() && "anchors must be loaded");

    SetBitPositions set_bits = get_row(row);
    SetBitPositions::iterator v = std::lower_bound(set_bits.begin(), set_bits.end(), column);
    return v != set_bits.end() && *v == column;
}


/**
 * Returns the given column.
 */
template <class BaseMatrix>
std::vector<BinaryMatrix::Row> RowDiff<BaseMatrix>::get_column(Column column) const {
    assert(anchor_.size() == diffs_.num_rows() && "anchors must be loaded");

    std::vector<Row> result;
    for (Row row = 0; row < num_rows(); ++row) {
        if (get(row, column))
            result.push_back(row);
    }
    return result;
}

template <class BaseMatrix>
BinaryMatrix::SetBitPositions RowDiff<BaseMatrix>::get_row(Row row) const {
    assert(anchor_.size() == diffs_.num_rows() && "anchors must be loaded");

    Vector<uint64_t> result = get_diff(row);
    std::sort(result.begin(), result.end());

    uint64_t boss_edge = graph_->kmer_to_boss_index(
            graph::AnnotatedSequenceGraph::anno_to_graph_index(row));
    const graph::boss::BOSS &boss = graph_->get_boss();

    while (!anchor_[row]) {
        graph::boss::BOSS::TAlphabet w = boss.get_W(boss_edge);
        assert(boss_edge > 1 && w != 0);

        // fwd always selects the last outgoing edge for a given node
        boss_edge = boss.fwd(boss_edge, w % boss.alph_size);
        row = graph::AnnotatedSequenceGraph::graph_to_anno_index(
                graph_->boss_to_kmer_index(boss_edge));

        auto diff_row = get_diff(row);
        std::sort(diff_row.begin(), diff_row.end());
        merge(&result, diff_row);
    }

    return result;
}

template <class BaseMatrix>
inline bool RowDiff<BaseMatrix>::load(std::istream &f) {
    if constexpr (!std::is_same_v<BaseMatrix, ColumnMajor>) {
        anchor_.load(f);
    }
    return diffs_.load(f);
}

template <class BaseMatrix>
inline void RowDiff<BaseMatrix>::serialize(std::ostream &f) const {
    if constexpr (!std::is_same_v<BaseMatrix, ColumnMajor>) {
        anchor_.serialize(f);
    }
    diffs_.serialize(f);
}

template <class BaseMatrix>
void RowDiff<BaseMatrix>::merge(Vector<uint64_t> *result, const Vector<uint64_t> &diff2) {
    assert(std::is_sorted(result->begin(), result->end()));
    assert(std::is_sorted(diff2.begin(), diff2.end()));

    if (diff2.empty())
        return;

    Vector<uint64_t> diff1;
    std::swap(*result, diff1);
    result->reserve(std::max(diff1.size(), diff2.size()));
    uint64_t idx1 = 0;
    uint64_t idx2 = 0;
    while (idx1 < diff1.size() && idx2 < diff2.size()) {
        if (diff1[idx1] == diff2[idx2]) {
            idx1++;
            idx2++;
        } else if (diff1[idx1] < diff2[idx2]) {
            result->push_back(diff1[idx1++]);
        } else {
            result->push_back(diff2[idx2++]);
        }
    }
    while (idx1 < diff1.size()) {
        result->push_back(diff1[idx1++]);
    }
    while (idx2 < diff2.size()) {
        result->push_back(diff2[idx2++]);
    }

    assert(std::is_sorted(result->begin(), result->end()));
}

} // namespace binmat
} // namespace annot
} // namespace mtg
