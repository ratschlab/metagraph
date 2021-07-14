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
    ColumnMajor(std::vector<std::unique_ptr<bit_vector>>&& columns)
        : columns_(std::move(columns)) {}

    uint64_t num_columns() const override { return columns_.size(); }
    uint64_t num_rows() const override;

    bool get(Row row, Column column) const override;
    SetBitPositions get_row(Row row) const override;
    std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const override;
    // query row and get ranks of each set bit in its column
    Vector<std::pair<Column, uint64_t>> get_column_ranks(Row row) const;
    std::vector<Vector<std::pair<Column, uint64_t>>>
    get_column_ranks(const std::vector<Row> &rows) const;
    std::vector<Row> get_column(Column column) const override;
    // get all selected rows appended with -1 and concatenated
    std::vector<Column> slice_rows(const std::vector<Row> &rows) const override;

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    // number of ones in the matrix
    uint64_t num_relations() const override;

    // Return all columns for which counts are greater than or equal to |min_count|.
    // Stop counting if count is greater than |count_cap|.
    std::vector<std::pair<Column, size_t /* count */>>
    sum_rows(const std::vector<std::pair<Row, size_t>> &index_counts,
             size_t min_count = 1,
             size_t count_cap = std::numeric_limits<size_t>::max()) const override;

    auto& data() { return columns_; }
    const auto& data() const { return columns_; }

  private:
    std::vector<std::unique_ptr<bit_vector>> columns_;
};

} // namespace binmat
} // namespace annot
} // namespace mtg

#endif // __COLUMN_MAJOR_HPP__
