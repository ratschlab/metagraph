#pragma once

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/vector_map.hpp"
#include "common/vector_set.hpp"
#include "common/vector.hpp"
#include "common/logger.hpp"
#include "common/utils/template_utils.hpp"
#include "common/hashers/hash.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace annot {

graph::DeBruijnGraph::node_index row_diff_successor(const graph::DeBruijnGraph &graph,
                                                    graph::DeBruijnGraph::node_index node,
                                                    const bit_vector &rd_succ);

namespace matrix {

const std::string kRowDiffAnchorExt = ".anchors";
const std::string kRowDiffForkSuccExt = ".rd_succ";

const size_t RD_PATH_RESERVE_SIZE = 2;


class IRowDiff {
  public:
    typedef bit_vector_small anchor_bv_type;
    typedef bit_vector_small fork_succ_bv_type;

    virtual ~IRowDiff() {}

    const graph::DeBruijnGraph* graph() const { return graph_; }
    void set_graph(const graph::DeBruijnGraph *graph) { graph_ = graph; }

    void load_fork_succ(const std::string &filename);
    void load_anchor(const std::string &filename);

    const anchor_bv_type& anchor() const { return anchor_; }

    const fork_succ_bv_type& fork_succ() const { return fork_succ_; }

  protected:
    // get row-diff paths starting at |row_ids|
    std::tuple<std::vector<BinaryMatrix::Row>, std::vector<std::vector<size_t>>, std::vector<size_t>>
    get_rd_ids(const std::vector<BinaryMatrix::Row> &row_ids, size_t num_threads = 1) const;

    template <class F, class G, class H, class Callback>
    void call_rows(const std::vector<BinaryMatrix::Row> &row_ids, F call_rd_rows, G add_diff,
                   H decode_diffs, Callback call_row, size_t num_threads) const;

    const graph::DeBruijnGraph *graph_ = nullptr;
    anchor_bv_type anchor_;
    fork_succ_bv_type fork_succ_;
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
    template <typename... Args>
    RowDiff(const graph::DeBruijnGraph *graph = nullptr, Args&&... args)
        : diffs_(std::forward<Args>(args)...) { graph_ = graph; }

    /**
     * Returns the number of set bits in the row-diff transformed matrix.
     */
    uint64_t num_relations() const override { return diffs_.num_relations(); }
    uint64_t num_columns() const override { return diffs_.num_columns(); }
    uint64_t num_rows() const override { return diffs_.num_rows(); }

    /**
     * Returns the given column.
     */
    std::vector<Row> get_column(Column column) const override;

    std::vector<SetBitPositions> get_rows(const std::vector<Row> &row_ids) const override;

    std::vector<SetBitPositions>
    get_rows_dict(std::vector<Row> *rows, size_t num_threads) const override;

    bool load(std::istream &f) override;
    void serialize(std::ostream &f) const override;

    const BaseMatrix& diffs() const { return diffs_; }
    BaseMatrix& diffs() { return diffs_; }

  private:
    static void add_diff(const SetBitPositions &diff, SetBitPositions *row);

    BaseMatrix diffs_;
};


/**
 * Returns the given column.
 */
template <class BaseMatrix>
std::vector<BinaryMatrix::Row> RowDiff<BaseMatrix>::get_column(Column column) const {
    assert(graph_ && "graph must be loaded");
    assert(anchor_.size() == diffs_.num_rows() && "anchors must be loaded");

    assert(!fork_succ_.size() || fork_succ_.size() == graph_->max_index() + 1);

    std::vector<Row> result;
    // TODO: implement a more efficient algorithm
    graph_->call_nodes([&](auto node) {
        auto row = graph::AnnotatedDBG::graph_to_anno_index(node);
        SetBitPositions set_bits = get_rows({ row })[0];
        if (std::binary_search(set_bits.begin(), set_bits.end(), column))
            result.push_back(row);
    });
    return result;
}

template <class F, class G, class H, class Callback>
void IRowDiff::call_rows(const std::vector<BinaryMatrix::Row> &row_ids, F call_rd_rows,
                         G add_diff, H decode_diffs, Callback call_row, size_t num_threads) const {
    assert(graph_ && "graph must be loaded");
    assert(anchor_.size() == graph::AnnotatedDBG::graph_to_anno_index(graph_->num_nodes() + 1)
                && "anchors must be loaded");
    assert(!fork_succ_.size() || fork_succ_.size() == graph_->max_index() + 1);

    // No sorting in order not to break the topological order for row-diff annotation

    // get row-diff paths
    auto [rd_ids, rd_paths_trunc, times_traversed] = get_rd_ids(row_ids, num_threads);

    auto rd_rows = call_rd_rows(rd_ids, num_threads);
    DEBUG_LOG("Queried batch of {} diffed rows", rd_ids.size());
    rd_ids = {};

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 1000)
    for (size_t i = 0; i < rd_rows.size(); ++i) {
        decode_diffs(&rd_rows[i]);
        std::sort(rd_rows[i].begin(), rd_rows[i].end());
    }

    // reconstruct annotation rows from row-diff
    typename decltype(rd_rows)::value_type result;

    for (size_t i = 0; i < row_ids.size(); ++i) {
        result.resize(0);
        // propagate back and reconstruct full annotations for predecessors
        for (auto it = rd_paths_trunc[i].rbegin(); it != rd_paths_trunc[i].rend(); ++it) {
            add_diff(rd_rows[*it], &result);
            // replace diff row with full reconstructed annotation
            if (--times_traversed[*it]) {
                rd_rows[*it] = result;
            } else {
                // free memory
                rd_rows[*it] = {};
            }
        }
        call_row(result);
    }
    DEBUG_LOG("Reconstructed annotations for {} rows", row_ids.size());
    assert(times_traversed == std::vector<size_t>(rd_rows.size(), 0));
}

template <class BaseMatrix>
std::vector<BinaryMatrix::SetBitPositions>
RowDiff<BaseMatrix>::get_rows(const std::vector<Row> &row_ids) const {
    std::vector<SetBitPositions> rows;
    call_rows(row_ids,
        [this](const std::vector<Row> &rd_ids, size_t num_threads) {
            return diffs_.get_rows(rd_ids, num_threads);
        },
        add_diff, [](SetBitPositions *row) {},
        [&](const SetBitPositions &row) { rows.push_back(row); },
        1
    );
    return rows;
}

template <class BaseMatrix>
std::vector<BinaryMatrix::SetBitPositions>
RowDiff<BaseMatrix>::get_rows_dict(std::vector<Row> *rows, size_t num_threads) const {
    VectorSet<SetBitPositions, utils::VectorHash> unique_rows;
    size_t i = 0;
    call_rows(*rows,
        [this](const std::vector<Row> &rd_ids, size_t num_threads) {
            return diffs_.get_rows(rd_ids, num_threads);
        },
        add_diff, [](SetBitPositions *row) {},
        [&](const SetBitPositions &row) {
            auto it = unique_rows.emplace(row).first;
            (*rows)[i++] = it - unique_rows.begin();
        },
        num_threads
    );
    return to_vector(std::move(unique_rows));
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
void RowDiff<BaseMatrix>::add_diff(const SetBitPositions &diff, SetBitPositions *row) {
    assert(std::is_sorted(row->begin(), row->end()));
    assert(std::is_sorted(diff.begin(), diff.end()));

    if (diff.empty())
        return;

    SetBitPositions result;
    result.reserve(row->size() + diff.size());
    std::set_symmetric_difference(row->begin(), row->end(),
                                  diff.begin(), diff.end(),
                                  std::back_inserter(result));
    row->swap(result);
}

} // namespace matrix
} // namespace annot
} // namespace mtg
