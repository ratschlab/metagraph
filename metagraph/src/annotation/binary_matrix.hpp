#ifndef __SPARSE_MATRIX_HPP__
#define __SPARSE_MATRIX_HPP__

#include <vector>

#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>

#include "utils/vectors.hpp"


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

    virtual bool load(std::istream &in) = 0;
    virtual void serialize(std::ostream &out) const = 0;

    // number of ones in the matrix
    virtual uint64_t num_relations() const = 0;
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

#endif // __SPARSE_MATRIX_HPP__
