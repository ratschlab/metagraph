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

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    // number of ones in the matrix
    uint64_t num_relations() const override;

    void point_to(const std::vector<std::unique_ptr<bit_vector>> &columns) {
        columns_ = &columns;
    }

    const auto& data() const { return *columns_; }

    using ColumnCallback = std::function<void(uint64_t, std::unique_ptr<bit_vector> &&)>;
    /**
     * Calls #callback on each of the columns. Note that this function "consumes" the
     * columns and leaves the object empty.
     */
    void call_columns(const ColumnCallback &callback);

  private:
    std::vector<std::unique_ptr<bit_vector>> data_;
    const std::vector<std::unique_ptr<bit_vector>> *columns_ = &data_;
};

} // namespace binmat
} // namespace annot
} // namespace mtg

#endif // __COLUMN_MAJOR_HPP__
