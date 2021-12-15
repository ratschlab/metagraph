#ifndef __CSR_MATRIX_HPP__
#define __CSR_MATRIX_HPP__

#include <vector>

#include "annotation/int_matrix/base/int_matrix.hpp"


namespace mtg {
namespace annot {
namespace matrix {

/**
 * Compressed Sparse Row Matrix
 *
 * Matrix which stores the non-zero values in row-major order.
 */
class CSRMatrix : public IntMatrix {
  public:
    explicit CSRMatrix(uint64_t num_rows = 0) : vector_(num_rows) {}

    CSRMatrix(Vector<RowValues>&& rows, uint64_t num_columns);

    // row is in [0, num_rows), column is in [0, num_columns)
    RowValues get_row_values(Row row) const { return vector_[row]; }

    std::vector<RowValues>
    get_row_values(const std::vector<Row> &rows) const;

    uint64_t num_columns() const { return num_columns_; }
    uint64_t num_rows() const { return vector_.size(); }
    uint64_t num_relations() const;

    // row is in [0, num_rows), column is in [0, num_columns)
    bool get(Row row, Column column) const;
    SetBitPositions get_row(Row row) const;
    std::vector<Row> get_column(Column column) const;

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

  private:
    uint64_t num_columns_ = 0;
    Vector<RowValues> vector_;
};

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __CSR_MATRIX_HPP__
