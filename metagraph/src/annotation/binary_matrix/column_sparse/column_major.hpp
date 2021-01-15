#ifndef __COLUMN_MAJOR_HPP__
#define __COLUMN_MAJOR_HPP__

#include <memory>

#include "common/vectors/bit_vector.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"


namespace mtg {
namespace annot {
namespace binmat {

class ColumnMajor : public BinaryMatrix {
  public:
    ColumnMajor() {}
    ColumnMajor(std::vector<std::unique_ptr<bit_vector>>&& columns);
    ColumnMajor(ColumnMajor&& other);

    uint64_t num_columns() const override { return columns_->size(); }
    uint64_t num_rows() const override;

    bool get(Row row, Column column) const override;
    SetBitPositions get_row(Row row) const override;
    std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const override;
    std::vector<Row> get_column(Column column) const override;
    // get all selected rows appended with -1 and concatenated
    std::vector<Column> slice_rows(const std::vector<Row> &rows) const override;

    void slice_columns(const std::vector<Column> &columns,
                       const std::function<void(Column, bitmap&&)> &callback,
                       size_t num_threads = 1) const override;

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    // number of ones in the matrix
    uint64_t num_relations() const override;

    void update_pointer(const std::vector<std::unique_ptr<bit_vector>> &columns) {
        columns_ = &columns;
    }

    const auto& data() const { return *columns_; }

    /**
     * Returns the columns vector and pilfers the existing columns.
     */
    std::vector<std::unique_ptr<bit_vector>> release_columns() {
        return std::move(data_);
    }

  private:
    std::vector<std::unique_ptr<bit_vector>> data_;
    const std::vector<std::unique_ptr<bit_vector>> *columns_ = &data_;
};

} // namespace binmat
} // namespace annot
} // namespace mtg

#endif // __COLUMN_MAJOR_HPP__
