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
#include "common/unix_tools.hpp"
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
    // Returns: (rd_ids, rd_paths_trunc, times_traversed, groups)
    // groups[g] records indices of row_ids paths traced in group g.
    // Rows from different groups access disjoint rd_rows entries, enabling
    // parallel reconstruction grouped by thread.
    std::tuple<std::vector<BinaryMatrix::Row>,
               std::vector<std::vector<size_t>>,
               std::vector<size_t>,
               std::vector<std::vector<size_t>>>
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
    using base_matrix_type = BaseMatrix;

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

    /**
     * Return rows (in arbitrary order) and update the row indexes in |rows|
     * to point to their respective rows in the vector returned.
     *
     * In contrast to get_rows_dict in most other classes, the rows here
     * are not deduplicated. Benchmarks on real queries showed that merging
     * identical reconstructed rows rarely shrinks the batch size by much
     * (only a few percent on dense annotation rows; occasionally on the order
     * of ~40% when duplication was high), while hashing requires a critical
     * section that dominates the query time. Hence, we skip deduplication.
     */
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
    assert(diffs_.num_rows() == graph_->max_index());
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
    assert(anchor_.size() == graph_->max_index() && "anchors must be loaded");
    assert(!fork_succ_.size() || fork_succ_.size() == graph_->max_index() + 1);

    if (row_ids.empty())
        return;

    // No sorting in order not to break the topological order for row-diff annotation

    // get row-diff paths
    Timer timer;
    // Unique row-diff row IDs to fetch and decode
    std::vector<BinaryMatrix::Row> rd_ids;
    // Truncated reconstruction paths (indices into rd_ids) per queried row
    std::vector<std::vector<size_t>> rd_paths_trunc;
    // Multiplicity for each queried row path returned by get_rd_ids()
    std::vector<size_t> times_traversed;
    // Independent groups of row paths that can be reconstructed in parallel
    std::vector<std::vector<size_t>> groups;
    std::tie(rd_ids, rd_paths_trunc, times_traversed, groups) = get_rd_ids(row_ids, num_threads);
    double rd_traversal_time = timer.elapsed();
    timer.reset();

    auto rd_rows = call_rd_rows(rd_ids, num_threads);
    double call_rd_rows_time = timer.elapsed();
    timer.reset();
    std::vector<BinaryMatrix::Row>().swap(rd_ids);

    size_t num_rd_bits = 0;
    size_t total_capacity = 0;
    // 200 rows per task, to make the task dispatch time negligible
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 200) \
        reduction(+:num_rd_bits,total_capacity)
    for (size_t i = 0; i < rd_rows.size(); ++i) {
        decode_diffs(&rd_rows[i]);
        std::sort(rd_rows[i].begin(), rd_rows[i].end(), utils::LessFirst());
        num_rd_bits += rd_rows[i].size();
        total_capacity += rd_rows[i].capacity();
    }

    double decode_diffs_time = timer.elapsed();
    timer.reset();

    // Reconstruct annotation rows from row-diff.
    // Since get_rd_ids skips cross-thread deduplication, rows from different
    // groups access disjoint rd_rows entries. We use that to reconstruct
    // each group in parallel, while preserving the sequential dependency
    // order within each group.
    using RowType = typename std::decay_t<decltype(rd_rows)>::value_type;

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (size_t g = 0; g < groups.size(); ++g) {
        RowType result;
        for (size_t i : groups[g]) {
            auto it = rd_paths_trunc[i].rbegin();
            result = rd_rows[*it];
            if (!(--times_traversed[*it]))
                RowType().swap(rd_rows[*it]);
            for (++it ; it != rd_paths_trunc[i].rend(); ++it) {
                add_diff(rd_rows[*it], &result);
                if (--times_traversed[*it]) {
                    rd_rows[*it] = result;
                } else {
                    RowType().swap(rd_rows[*it]);
                }
            }
            call_row(i, result);
        }
    }

    common::logger->trace("RD query [threads: {}, rows: {} -> {} ({:.1f}x)] -- "
            "traversal: {:.2f} sec, call_rd_rows: {:.2f} sec (set bits: {}, capacity: {}), "
            "decoding: {:.2f} sec, reconstruction: {:.2f} sec",
            num_threads, row_ids.size(), rd_rows.size(), (double)rd_rows.size()/row_ids.size(),
            rd_traversal_time, call_rd_rows_time, num_rd_bits, total_capacity,
            decode_diffs_time, timer.elapsed());

    assert(times_traversed == std::vector<size_t>(rd_rows.size(), 0));
    assert(std::all_of(rd_rows.begin(), rd_rows.end(), [](const auto &v) { return v.empty(); }));
}

template <class BaseMatrix>
std::vector<BinaryMatrix::SetBitPositions>
RowDiff<BaseMatrix>::get_rows(const std::vector<Row> &row_ids) const {
    std::vector<SetBitPositions> rows(row_ids.size());
    call_rows(row_ids,
        [this](const std::vector<Row> &rd_ids, size_t num_threads) {
            return diffs_.get_rows(rd_ids, num_threads);
        },
        add_diff, [](SetBitPositions *row) {},
        [&](size_t i, const SetBitPositions &row) { rows[i] = row; },
        1
    );
    return rows;
}

template <class BaseMatrix>
std::vector<BinaryMatrix::SetBitPositions>
RowDiff<BaseMatrix>::get_rows_dict(std::vector<Row> *rows, size_t num_threads) const {
    std::vector<SetBitPositions> rows_dict(rows->size());
    call_rows(*rows,
        [this](const std::vector<Row> &rd_ids, size_t num_threads) {
            return diffs_.get_rows(rd_ids, num_threads);
        },
        add_diff, [](SetBitPositions *row) {},
        [&](size_t i, const SetBitPositions &row) {
            rows_dict[i] = row;
            (*rows)[i] = i;
        },
        num_threads
    );
    return rows_dict;
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
