#ifndef __BRWT_HPP__
#define __BRWT_HPP__

#include <vector>
#include <unordered_map>
#include <memory>

#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/range_partition.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"


namespace mtg {
namespace annot {
namespace matrix {

// The Multi-BRWT compressed binary matrix representation
class BRWT : public BinaryMatrix, public GetEntrySupport {
    friend class BRWTBuilder;
    friend class BRWTBottomUpBuilder;
    friend class BRWTOptimizer;

    typedef uint32_t Child;

  public:
    BRWT() : nonzero_rows_(new bit_vector_smallrank()) {}

    uint64_t num_columns() const override { return assignments_.size(); }
    uint64_t num_rows() const override { return nonzero_rows_->size(); }

    bool get(Row row, Column column) const override;
    std::vector<Row> get_column(Column column) const override;
    std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const override;
    std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows, size_t num_threads) const override;
    // query row and get ranks of each set bit in its column
    std::vector<Vector<std::pair<Column, uint64_t>>>
    get_column_ranks(const std::vector<Row> &rows, size_t num_threads) const;
    void call_rows(const std::function<void(const SetBitPositions &)> &callback,
                   bool show_progress = common::get_verbose()) const;

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    // number of ones in the matrix
    uint64_t num_relations() const override;

    // internal stats functions
    double avg_arity() const;
    uint64_t num_nodes() const;
    double shrinking_rate() const;

    void print_tree_structure(std::ostream &os) const;

  private:
    // breadth-first traversal
    void BFT(std::function<void(const BRWT &node)> callback) const;
    // get all selected rows appended with -1 and concatenated
    // helper function for querying rows in batches
    // appends to `slice`
    template <typename T>
    void slice_rows(const std::vector<Row> &rows, Vector<T> *slice,
                    ThreadPool *thread_pool = nullptr) const;

    template <typename RowT>
    std::vector<RowT>
    slice_rows_parallel(const std::vector<Row> &row_ids, size_t num_threads) const;

    // `call_stack` tracks the child indexes at each level during parallel tree
    // traversal, so that column indices can be remapped when merging results from subtrees.
    template <typename T>
    void slice_rows(const std::vector<Row> &row_ids, std::vector<size_t> rows,
                    const BRWT *root, std::vector<size_t> call_stack,
                    size_t max_columns_cutoff, ThreadPool &thread_pool,
                    std::function<void(std::vector<size_t>&&, Vector<T>&&)> call_slice) const;

    // Returns (`nonzero_indices`, `child_row_ids`): indices into `rows` for
    // the rows with the set bit in `nonzero_rows_`, and their recalculated
    // (via rank) row indices for the respective rows in the children.
    std::pair<std::vector<size_t>, std::vector<BRWT::Row>>
    get_nonzero_rows(const std::vector<Row> &rows,
                     ThreadPool *thread_pool = nullptr,
                     bool adaptive_chunk_size = true) const;

    // assigns columns to the child nodes
    RangePartition assignments_;
    std::unique_ptr<bit_vector> nonzero_rows_;
    // generally, these child matrices can be abstract BinaryMatrix instances
    std::vector<std::unique_ptr<BRWT>> child_nodes_;
};

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __BRWT_HPP__
