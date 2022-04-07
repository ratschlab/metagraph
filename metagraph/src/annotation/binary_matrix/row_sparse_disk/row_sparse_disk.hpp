#pragma once

#include <iostream>
#include <string>
#include <vector>

//#include <sdsl/vlc_vector.hpp>
//#include <sdsl/vlc_vector_buffer.hpp>

#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"

namespace mtg {
namespace annot {
namespace binmat {


class RowSparseDisk : public BinaryMatrix {
  public:
    RowSparseDisk() {}
    RowSparseDisk(const std::function<void(const RowCallback&)> &call_rows,
              uint64_t num_columns,
              uint64_t num_rows,
              uint64_t num_relations);

    uint64_t num_columns() const override { return num_columns_; }
    uint64_t num_rows() const override { return num_rows_; }

    bool get(Row row, Column column) const override;
    SetBitPositions get_row(Row row) const override;
    std::vector<Row> get_column(Column column) const override;

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    // number of ones in the matrix
    uint64_t num_relations() const override { return set_bits_.size(); }

    static void serialize(const std::function<void(binmat::BinaryMatrix::RowCallback)> &call_rows,
                          const std::string& filename, uint64_t num_cols, uint64_t num_set_bits, uint64_t num_rows);
  private:
    //sdsl::vlc_vector<> set_bits_;
    sdsl::int_vector_buffer<> set_bits_;
    bit_vector_small boundary_;
    uint64_t num_columns_ = 0;
    uint64_t num_rows_ = 0;
};


} // namespace binmat
} // namespace annot
} // namespace mtg
