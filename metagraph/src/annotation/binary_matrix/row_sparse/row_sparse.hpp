#pragma once

#include <iostream>
#include <string>
#include <vector>

#include <sdsl/packed_vector_view.hpp>
#include <sdsl/vlc_vector.hpp>

#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "annotation/binary_matrix/base/decode_budget.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"

namespace mtg {
namespace annot {
namespace matrix {

class RowSparse : public RowMajor {
  public:
    RowSparse() {}
    RowSparse(const std::function<void(const RowCallback&)> &call_rows,
              uint64_t num_columns,
              uint64_t num_rows,
              uint64_t num_relations);

    uint64_t num_columns() const override { return num_columns_; }
    uint64_t num_rows() const override { return num_rows_; }

    SetBitPositions get_row(Row row) const override;
    std::vector<Row> get_column(Column column) const override;
    // The budget-aware decode path (decode_budget.hpp, as BRWT::row_columns): the set bits
    // of one row, ascending, appended to |scratch|. The row is one run of the delta-coded
    // set bits, so one read; every growth of the scratch is charged before it happens.
    // false: a growth did not fit (|scratch| keeps what was appended before).
    bool row_columns(Row row, ChargedBuffer<Column> *scratch, DecodeBudget &budget) const;

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    // number of ones in the matrix
    uint64_t num_relations() const override { return set_bits_.size(); }

  private:
    sdsl::vlc_vector<sdsl::coder::elias_delta<>, 128, 0, sdsl::packed_vector_view> set_bits_;
    bit_vector_small boundary_;
    uint64_t num_columns_ = 0;
    uint64_t num_rows_ = 0;
};

} // namespace matrix
} // namespace annot
} // namespace mtg
