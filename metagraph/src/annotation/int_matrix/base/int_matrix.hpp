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
    virtual RowValues
    sum_row_values(const std::vector<std::pair<Row, size_t>> &index_counts,
                   size_t min = 1,
                   size_t cap = std::numeric_limits<size_t>::max()) const;
};


// Entries are tuples and aggregated values are tuple sizes
class MultiIntMatrix : public IntMatrix {
  public:
    typedef SmallVector<uint64_t> Tuple;
    typedef Vector<std::pair<Column, Tuple>> RowTuples;

    virtual ~MultiIntMatrix() {}

    // return tuple sizes (if not zero) at each entry
    virtual RowValues get_row_values(Row row) const;

    virtual std::vector<RowValues>
    get_row_values(const std::vector<Row> &rows) const;

    // return total number of attributes in all tuples
    virtual uint64_t num_attributes() const = 0;

    // return entries of the matrix -- where each entry is a set of integers
    virtual RowTuples get_row_tuples(Row row) const = 0;

    virtual std::vector<RowTuples>
    get_row_tuples(const std::vector<Row> &rows) const = 0;

    virtual bool load_tuples(std::istream &in) = 0;
    virtual void serialize_tuples(std::ostream &out) const = 0;
};

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __INT_MATRIX_HPP__
