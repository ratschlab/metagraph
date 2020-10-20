#pragma once

#include <fstream>
#include <string>
#include <vector>

#include <sdsl/enc_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>
#include <sdsl/rrr_vector.hpp>

#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "common/logger.hpp"
#include "common/vector.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"

namespace mtg {
namespace annot {
namespace binmat {

class RowSparse : public BinaryMatrix {
  public:
    RowSparse() {}
    RowSparse(const std::function<void(const RowCallback&)> &call_rows,
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
    uint64_t num_relations() const override { return elements_.size(); }

  private:
    sdsl::enc_vector<> elements_;
    sdsl::rrr_vector<> boundary_;
    sdsl::rrr_vector<>::select_1_type sboundary_;
    uint64_t num_columns_ = 0;
    uint64_t num_rows_ = 0;
};

} // namespace binmat
} // namespace annot
} // namespace mtg
