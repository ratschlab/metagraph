#ifndef __INT_MATRIX_HPP__
#define __INT_MATRIX_HPP__

#include <vector>

#include <sdsl/int_vector.hpp>

#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "common/vector_map.hpp"


namespace mtg {
namespace annot {
namespace matrix {

class IntMatrix {
  public:
    typedef uint64_t Value;
    typedef Vector<std::pair<BinaryMatrix::Column, Value>> RowValues;

    virtual ~IntMatrix() {}

    // |row| is in [0, num_rows), |column| is in [0, num_columns)
    virtual std::vector<RowValues> get_row_values(const std::vector<BinaryMatrix::Row> &rows) const = 0;

    virtual void call_row_values(const std::function<void(uint64_t, RowValues&&, size_t)> &callback,
                                 bool ordered = true) const;

    virtual std::vector<VectorMap<uint64_t, size_t>> get_histograms(const std::vector<size_t> &min_counts = {},
                                                                    sdsl::bit_vector *unmark_discarded = nullptr) const;

    // sum up values for each column with at least |min_count| non-zero values
    virtual RowValues
    sum_row_values(const std::vector<std::pair<BinaryMatrix::Row, size_t>> &index_counts,
                   size_t min_count = 1) const;

    virtual const BinaryMatrix& get_binary_matrix() const = 0;
};


// Entries are tuples and their aggregated `values` are tuple sizes
class MultiIntMatrix : public IntMatrix {
  public:
    typedef SmallVector<uint64_t> Tuple;
    typedef Vector<std::pair<BinaryMatrix::Column, Tuple>> RowTuples;

    virtual ~MultiIntMatrix() {}

    // return tuple sizes (if not zero) at each entry
    virtual std::vector<RowValues> get_row_values(const std::vector<BinaryMatrix::Row> &rows) const;

    // return total number of attributes in all tuples
    virtual uint64_t num_attributes() const = 0;

    // return entries of the matrix -- where each entry is a set of integers
    virtual std::vector<RowTuples> get_row_tuples(const std::vector<BinaryMatrix::Row> &rows) const = 0;
};

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __INT_MATRIX_HPP__
