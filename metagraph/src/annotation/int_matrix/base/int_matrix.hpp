#ifndef __INT_MATRIX_HPP__
#define __INT_MATRIX_HPP__

#include <vector>

#include "annotation/binary_matrix/base/binary_matrix.hpp"


namespace mtg {
namespace annot {
namespace matrix {

class IntMatrix : public binmat::BinaryMatrix {
  public:
    typedef uint64_t Value;
    typedef Vector<std::pair<Column, Value>> RowValues;

    virtual ~IntMatrix() {}

    // |row| is in [0, num_rows), |column| is in [0, num_columns)
    virtual RowValues get_row_values(Row row) const = 0;

    virtual std::vector<RowValues>
    get_row_values(const std::vector<Row> &rows) const = 0;

    // Get all columns for which the sum of the values in queried rows
    // is greater than or equal to |min|. Stop counting if the sum is
    // greater than |cap|.
    virtual std::vector<std::pair<Column, size_t /* count */>>
    sum_row_values(const std::vector<std::pair<Row, size_t>> &index_counts,
                   size_t min = 1,
                   size_t cap = std::numeric_limits<size_t>::max()) const;
};

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __INT_MATRIX_HPP__
