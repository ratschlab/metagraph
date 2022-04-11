#ifndef __BRWT_HPP__
#define __BRWT_HPP__

#include <vector>
#include <unordered_map>
#include <memory>

#include <progress_bar.hpp>
#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/range_partition.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"


namespace mtg {
namespace annot {
namespace binmat {

// The Multi-BRWT compressed binary matrix representation
class BRWT : public BinaryMatrix {
    friend class BRWTBuilder;
    friend class BRWTBottomUpBuilder;
    friend class BRWTOptimizer;

    typedef uint32_t Child;

    void build_histo(std::vector<size_t>& cnts_per_row, const std::vector<uint64_t>& row_ids_mapping, ProgressBar& progress_bar) const;

  public:
    BRWT() : nonzero_rows_(new bit_vector_smallrank()) {}

    uint64_t num_columns() const override { return assignments_.size(); }
    uint64_t num_rows() const override { return nonzero_rows_->size(); }

    std::vector<size_t> get_rows_set_bits_histo() const;

    bool get(Row row, Column column) const override;
    SetBitPositions get_row(Row row) const override;
    std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const override;
    std::vector<Row> get_column(Column column) const override;
    // get all selected rows appended with -1 and concatenated
    std::vector<Column> slice_rows(const std::vector<Row> &rows) const override;
    // query row and get ranks of each set bit in its column
    Vector<std::pair<Column, uint64_t>> get_column_ranks(Row row) const;
    std::vector<Vector<std::pair<Column, uint64_t>>>
    get_column_ranks(const std::vector<Row> &rows) const;

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    // number of ones in the matrix
    uint64_t num_relations() const override;

    // internal stats functions
    double avg_arity() const;
    uint64_t num_nodes() const;
    double shrinking_rate() const;
    uint64_t total_column_size() const;
    uint64_t total_num_set_bits() const;

    void print_tree_structure(std::ostream &os) const;

  private:
    // breadth-first traversal
    void BFT(std::function<void(const BRWT &node)> callback) const;
    // helper function for querying rows in batches
    template <typename T>
    std::vector<T> slice_rows(const std::vector<Row> &rows) const;

    // assigns columns to the child nodes
    RangePartition assignments_;
    std::unique_ptr<bit_vector> nonzero_rows_;
    // generally, these child matrices can be abstract BinaryMatrix instances
    std::vector<std::unique_ptr<BRWT>> child_nodes_;
};

} // namespace binmat
} // namespace annot
} // namespace mtg

#endif // __BRWT_HPP__
