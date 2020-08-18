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
namespace binmat {

// The Multi-BRWT compressed binary matrix representation
class BRWT : public BinaryMatrix {
    friend class BRWTBuilder;
    friend class BRWTBottomUpBuilder;
    friend class BRWTOptimizer;

    typedef uint32_t Child;

  public:
    BRWT() : nonzero_rows_(new bit_vector_smallrank()) {}

    uint64_t num_columns() const { return assignments_.size(); }
    uint64_t num_rows() const { return nonzero_rows_->size(); }

    bool get(Row row, Column column) const;
    SetBitPositions get_row(Row row) const;
    std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const;
    std::vector<Row> get_column(Column column) const;

    // get all selected rows appended with -1 and concatenated
    std::vector<Column> slice_rows(const std::vector<Row> &rows) const;

    void slice_columns(const std::vector<Column> &columns,
                       const ColumnCallback &callback) const;

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    // number of ones in the matrix
    uint64_t num_relations() const;

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

    // assigns columns to the child nodes
    RangePartition assignments_;
    std::unique_ptr<bit_vector> nonzero_rows_;
    std::vector<std::unique_ptr<BinaryMatrix>> child_nodes_;
};

} // namespace binmat
} // namespace annot
} // namespace mtg

#endif // __BRWT_HPP__
