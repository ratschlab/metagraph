#ifndef __COLUMN_MAJOR_HPP__
#define __COLUMN_MAJOR_HPP__

#include <memory>

#include "common/vectors/bit_vector.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"


class ColumnMajor : public BinaryMatrix {
  public:
    ColumnMajor() {}
    ColumnMajor(std::vector<std::unique_ptr<bit_vector>>&& columns);

    uint64_t num_columns() const { return columns_->size(); }
    uint64_t num_rows() const;

    bool get(Row row, Column column) const;
    SetBitPositions get_row(Row row) const;
    std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const;
    std::vector<Row> get_column(Column column) const;
    // get all selected rows appended with -1 and concatenated
    std::vector<Column> slice_rows(const std::vector<Row> &rows) const;

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    // number of ones in the matrix
    uint64_t num_relations() const;

    static ColumnMajor
    construct_view(const std::vector<std::unique_ptr<bit_vector>> &columns) {
        ColumnMajor view;
        view.columns_ = &columns;
        return view;
    }

    const auto& data() const { return *columns_; }

  private:
    std::vector<std::unique_ptr<bit_vector>> data_;
    const std::vector<std::unique_ptr<bit_vector>> *columns_ = &data_;
};

#endif // __COLUMN_MAJOR_HPP__
