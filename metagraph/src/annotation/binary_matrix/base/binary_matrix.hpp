#ifndef __SPARSE_MATRIX_HPP__
#define __SPARSE_MATRIX_HPP__

#include <vector>

#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>

#include "common/vector.hpp"


class bitmap;

namespace mtg {
namespace annot {
namespace binmat {

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
    virtual bool get(Row row, Column column) const = 0;
    virtual SetBitPositions get_row(Row row) const = 0;
    virtual std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const;
    virtual std::vector<Row> get_column(Column column) const = 0;

    // get all selected rows appended with -1 and concatenated
    virtual std::vector<Column> slice_rows(const std::vector<Row> &rows) const;

    virtual void call_columns(const std::vector<Column> &columns,
                              const std::function<void(Column, const bitmap&)> &callback,
                              size_t num_threads = 1) const;

    virtual bool load(std::istream &in) = 0;
    virtual void serialize(std::ostream &out) const = 0;

    // number of ones in the matrix
    virtual uint64_t num_relations() const = 0;

    // Return all columns for which counts are greater than or equal to |min_count|.
    // Stop counting if count is greater than |count_cap|.
    virtual std::vector<std::pair<Column, size_t /* count */>>
    sum_rows(const std::vector<std::pair<Row, size_t>> &index_counts,
             size_t min_count = 1,
             size_t count_cap = std::numeric_limits<size_t>::max()) const;
};


class RainbowMatrix : public BinaryMatrix {
  public:
    virtual ~RainbowMatrix() {}

    using BinaryMatrix::get_rows;

    virtual SetBitPositions get_row(Row row) const final {
        return code_to_row(get_code(row));
    }

    // Return unique rows (in arbitrary order) and update the row indexes
    // in |rows| to point to their respective rows in the vector returned.
    virtual std::vector<SetBitPositions>
    get_rows(std::vector<Row> *rows, size_t num_threads = 1) const;

    // Return all columns for which counts are greater than or equal to |min_count|.
    // Stop counting if count is greater than |count_cap|.
    virtual std::vector<std::pair<Column, size_t /* count */>>
    sum_rows(const std::vector<std::pair<Row, size_t>> &index_counts,
             size_t min_count = 1,
             size_t count_cap = std::numeric_limits<size_t>::max()) const;

    virtual uint64_t num_distinct_rows() const = 0;

  private:
    virtual uint64_t get_code(Row row) const = 0;
    virtual SetBitPositions code_to_row(Row row) const = 0;
};


class BinaryMatrixRowDynamic : public BinaryMatrix {
  public:
    virtual ~BinaryMatrixRowDynamic() {}

    virtual void set(Row row, Column column) = 0;
    // fast set: do not check for existance before adding the label
    virtual void force_set(Row row, Column column) { set(row, column); }

    virtual void clear_row(Row row) = 0;
    virtual void insert_rows(const std::vector<Row> &rows) = 0;
};


// Row streamer -- read rows from a serialized row major binary matrix
template <typename RowType = BinaryMatrix::SetBitPositions>
class StreamRows {
  public:
    StreamRows(const std::string &filename, size_t offset);

    //TODO: implement constructor from stream once
    //      it's implemented for sdsl::int_vector_buffer<>.
    //      Then, use StreamRows to simplify load functions.
    // StreamRows(std::istream &instream);

    // return nullptr after all rows have been called
    RowType* next_row();

  private:
    RowType row_;
    sdsl::int_vector_buffer<> inbuf_;
    uint64_t i_ = 0;
};

// Write matrix to the end
void append_row_major(const std::string &filename,
                      const std::function<void(BinaryMatrix::RowCallback)> &call_rows,
                      uint64_t num_cols);

} // namespace binmat
} // namespace annot
} // namespace mtg

#endif // __SPARSE_MATRIX_HPP__
