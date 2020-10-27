#ifndef __BIN_REL_WT_HPP__
#define __BIN_REL_WT_HPP__

#include <brwt/binary_relation.h>

#include "annotation/binary_matrix/base/binary_matrix.hpp"

class bit_vector;


namespace mtg {
namespace annot {
namespace binmat {

class BinRelWT : public BinaryMatrix {
  public:
    BinRelWT() {}

    BinRelWT(const std::function<void(const RowCallback &)> &generate_rows,
             uint64_t num_set_bits, uint64_t num_columns);

    BinRelWT(std::vector<std::unique_ptr<bit_vector>>&& columns);

    uint64_t num_columns() const override;
    uint64_t num_rows() const override;

    bool get(Row row, Column column) const override;
    SetBitPositions get_row(Row row) const override;
    std::vector<Row> get_column(Column column) const override;

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    // number of ones in the matrix
    uint64_t num_relations() const override;

  private:
    bool is_zero_row(Row row) const;
    bool is_zero_column(Column column) const;

    brwt::binary_relation binary_relation_;
    uint64_t num_labels = 0;
    uint64_t max_used_label = 0;
    uint64_t max_used_object = 0;
    uint64_t num_objects = 0;
};

} // namespace binmat
} // namespace annot
} // namespace mtg

#endif // __BIN_REL_WT_HPP__
