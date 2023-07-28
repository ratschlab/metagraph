#ifndef __SPARSE_MATRIX_HPP__
#define __SPARSE_MATRIX_HPP__

#include <vector>
#include <functional>
#include <istream>
#include <ostream>

#include "common/vector.hpp"


class bitmap;

namespace mtg {
namespace annot {
namespace matrix {

class BinaryMatrix {
  public:
    typedef uint64_t Row;
    typedef uint64_t Column;

    typedef Vector<Column> SetBitPositions;
    typedef std::function<void(const SetBitPositions &)> RowCallback;
    typedef std::function<void(Row, Column)> ValueCallback;

    virtual ~BinaryMatrix() {}

    virtual uint64_t num_columns() const = 0;
    virtual uint64_t num_rows() const = 0;

    // row is in [0, num_rows), column is in [0, num_columns)
    virtual std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const = 0;
    // Return unique rows (in arbitrary order) and update the row indexes
    // in |rows| to point to their respective rows in the vector returned.
    virtual std::vector<SetBitPositions> get_rows_dict(std::vector<Row> *rows,
                                                       size_t num_threads = 1) const;
    virtual std::vector<Row> get_column(Column column) const = 0;

    // For each column id in columns, run callback on its respective index in columns
    // and a bitmap represnting the column
    virtual void call_columns(const std::vector<Column> &columns,
                              const std::function<void(size_t, const bitmap&)> &callback,
                              size_t num_threads = 1) const;

    virtual bool load(std::istream &in) = 0;
    virtual void serialize(std::ostream &out) const = 0;

    // number of ones in the matrix
    virtual uint64_t num_relations() const = 0;

    // Return all columns for which counts are greater than or equal to |min_count|.
    virtual std::vector<std::pair<Column, size_t /* count */>>
    sum_rows(const std::vector<std::pair<Row, size_t>> &index_counts,
             size_t min_count = 1) const;
};


class RainbowMatrix : public BinaryMatrix {
  public:
    virtual ~RainbowMatrix() {}

    // row is in [0, num_rows), column is in [0, num_columns)
    virtual std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const final;

    // Return unique rows (in arbitrary order) and update the row indexes
    // in |rows| to point to their respective rows in the vector returned.
    virtual std::vector<SetBitPositions> get_rows_dict(std::vector<Row> *rows,
                                                       size_t num_threads = 1) const final;

    // Return all columns for which counts are greater than or equal to |min_count|.
    virtual std::vector<std::pair<Column, size_t /* count */>>
    sum_rows(const std::vector<std::pair<Row, size_t>> &index_counts,
             size_t min_count = 1) const;

    virtual uint64_t num_distinct_rows() const = 0;

  private:
    virtual uint64_t get_code(Row row) const = 0;
    virtual std::vector<SetBitPositions> codes_to_rows(const std::vector<Row> &rows) const = 0;
};


class RowMajor : public BinaryMatrix {
  public:
    virtual ~RowMajor() {}

    virtual SetBitPositions get_row(Row row) const = 0;

    // row is in [0, num_rows), column is in [0, num_columns)
    virtual std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const final;
};

class GetEntrySupport {
  public:
    virtual ~GetEntrySupport() {}

    virtual bool get(BinaryMatrix::Row row, BinaryMatrix::Column column) const = 0;
};

class BinaryMatrixRowDynamic : public RowMajor {
  public:
    virtual ~BinaryMatrixRowDynamic() {}

    virtual void set(Row row, Column column) = 0;
    // fast set: do not check for existance before adding the label
    virtual void force_set(Row row, Column column) { set(row, column); }

    virtual void clear_row(Row row) = 0;
    virtual void insert_rows(const std::vector<Row> &rows) = 0;
};

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __SPARSE_MATRIX_HPP__
