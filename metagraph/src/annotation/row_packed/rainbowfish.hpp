#ifndef __RAINBOWFISH_HPP__
#define __RAINBOWFISH_HPP__

#include <vector>

#include "binary_matrix.hpp"
#include "bit_vector.hpp"
#include "flat_matrix.hpp"


class Rainbowfish : public BinaryMatrix {
  public:
    using RowCallback = std::function<void(const std::vector<Column> &)>;

    Rainbowfish() {}
    Rainbowfish(const std::function<void(RowCallback)> &call_rows,
                uint64_t num_columns,
                uint64_t buffer_size = static_cast<uint64_t>(-1));

    uint64_t num_columns() const { return num_columns_; }

    uint64_t num_rows() const;

    uint64_t num_distinct_rows() const;

    // row is in [0, num_rows), column is in [0, num_columns)
    bool get(Row row, Column column) const;
    std::vector<Column> get_row(Row row) const;
    std::vector<Row> get_column(Column column) const;

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    // number of ones in the matrix
    uint64_t num_relations() const;

    const bit_vector_rrr<>& get_row_codes() const { return row_codes_; }
    const bit_vector_rrr<>& get_row_code_delimiters() const { return row_code_delimiters_; }

    const std::vector<std::unique_ptr<BinaryMatrix>>& get_distinct_rows() const {
        return reduced_matrix_;
    }

  private:
    using ReducedMatrixType = RowConcatenated<bit_vector_small>;

    uint64_t num_columns_ = 0;
    uint64_t num_relations_ = 0;
    uint64_t buffer_size_ = static_cast<uint64_t>(-1);

    bit_vector_rrr<> row_codes_;
    bit_vector_rrr<> row_code_delimiters_;
    std::vector<std::unique_ptr<BinaryMatrix>> reduced_matrix_;

    uint64_t get_code(Row row) const;
};

#endif // __RAINBOWFISH_HPP__
