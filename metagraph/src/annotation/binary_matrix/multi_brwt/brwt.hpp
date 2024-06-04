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
    // query row and get ranks of each set bit in its column
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

    void print_tree_structure(std::ostream &os) const;

  private:
    // breadth-first traversal
    void BFT(std::function<void(const BRWT &node)> callback) const;
    // get all selected rows appended with -1 and concatenated
    // helper function for querying rows in batches
    // appends to `slice`
    template <typename T>
    void slice_rows(const std::vector<Row> &rows, Vector<T> *slice) const;

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
