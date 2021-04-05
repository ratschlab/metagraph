#ifndef __COUNT_MATRIX_HPP__
#define __COUNT_MATRIX_HPP__

#include <vector>

#include "annotation/binary_matrix/base/binary_matrix.hpp"


namespace mtg {
namespace annot {
namespace matrix {

class CountMatrix : public binmat::BinaryMatrix {
  public:
    typedef uint64_t Count;
    typedef Vector<std::pair<Column, Count>> ColumnCounts;

    virtual ~CountMatrix() {}

    // row is in [0, num_rows), column is in [0, num_columns)
    virtual ColumnCounts get_row_counts(Row row) const = 0;

    virtual std::vector<ColumnCounts>
    get_row_counts(const std::vector<Row> &rows) const = 0;

    // Return all columns for which sum counts are greater than or equal to
    // |min_count|. Stop counting if sum of counts is greater than |count_cap|.
    virtual std::vector<std::pair<Column, size_t /* count */>>
    sum_row_counts(const std::vector<std::pair<Row, size_t>> &index_counts,
                   size_t min_count = 1,
                   size_t count_cap = std::numeric_limits<size_t>::max()) const;
};

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __COUNT_MATRIX_HPP__
