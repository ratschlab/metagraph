#ifndef __COLUMN_MAJOR_HPP__
#define __COLUMN_MAJOR_HPP__

#include <memory>

#include "common/vectors/bit_vector.hpp"
#include "binary_matrix/base/binary_matrix.hpp"


class ColumnMajor : public BinaryMatrix {
  public:
    ColumnMajor() {}
    ColumnMajor(const std::vector<std::unique_ptr<bit_vector_sd>> &columns);
    ColumnMajor(std::vector<std::unique_ptr<bit_vector_sd>>&& columns);

    ColumnMajor(const ColumnMajor &other) = default;
    ColumnMajor& operator=(const ColumnMajor &other) = default;

    ColumnMajor(ColumnMajor&& other) = default;
    ColumnMajor& operator=(ColumnMajor&& other) = default;

    uint64_t num_columns() const { return columns_.size(); }
    uint64_t num_rows() const;

    bool get(Row row, Column column) const;
    SetBitPositions get_row(Row row) const;
    std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const;
    std::vector<Row> get_column(Column column) const;

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    // number of ones in the matrix
    uint64_t num_relations() const;

  private:
    std::vector<std::unique_ptr<bit_vector_sd>> columns_;
};

#endif // __COLUMN_MAJOR_HPP__
