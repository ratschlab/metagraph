#ifndef __VECTOR_ROW_BINMAT_HPP__
#define __VECTOR_ROW_BINMAT_HPP__

#include <vector>

#include "common/vector.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"


namespace mtg {
namespace annot {
namespace binmat {

template <typename RowType = SmallVector<uint32_t>>
class VectorRowBinMat : public BinaryMatrixRowDynamic {
  public:
    typedef RowType row_type;

    explicit VectorRowBinMat(uint64_t num_rows = 0) : vector_(num_rows) {}

    VectorRowBinMat(std::vector<RowType>&& rows, uint64_t num_columns);

    uint64_t num_columns() const { return num_columns_; }
    uint64_t num_rows() const { return vector_.size(); }

    bool get(Row row, Column column) const;
    SetBitPositions get_row(Row row) const;
    std::vector<Row> get_column(Column column) const;

    void set(Row row, Column column);
    // do not check for existance before adding the label
    void force_set(Row row, Column column);
    void clear_row(Row row);

    // sort all indexes in rows and leave only distinct ones
    void standardize_rows();

    void insert_rows(const std::vector<Row> &rows);

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    // number of ones in the matrix
    uint64_t num_relations() const;
    // matrix density
    double density() const;

  private:
    uint64_t num_columns_ = 0;
    std::vector<RowType> vector_;
};

} // namespace binmat
} // namespace anno
} // namespace mtg

#endif // __VECTOR_ROW_BINMAT_HPP__
