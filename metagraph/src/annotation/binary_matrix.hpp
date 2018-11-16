#ifndef __SPARSE_MATRIX_HPP__
#define __SPARSE_MATRIX_HPP__

#include <vector>

#include "bit_vector.hpp"


class BinaryMatrix {
  public:
    typedef uint64_t Row;
    typedef uint64_t Column;

    virtual ~BinaryMatrix() {}

    virtual uint64_t num_columns() const = 0;
    virtual uint64_t num_rows() const = 0;

    // row is in [0, num_rows), column is in [0, num_columns)
    virtual bool get(Row row, Column column) const = 0;
    virtual std::vector<Column> get_row(Row row) const = 0;
    virtual std::vector<Row> get_column(Column column) const = 0;

    virtual bool load(std::istream &in) = 0;
    virtual void serialize(std::ostream &out) const = 0;

    // number of ones in the matrix
    virtual uint64_t num_relations() const = 0;
};


class BinaryMatrixRowDynamic : public BinaryMatrix {
  public:
    virtual ~BinaryMatrixRowDynamic() {}

    virtual void set(Row row, Column column) = 0;
    virtual void clear_row(Row row) = 0;
    virtual void insert_rows(const std::vector<Row> &rows) = 0;
};


#endif // __SPARSE_MATRIX_HPP__