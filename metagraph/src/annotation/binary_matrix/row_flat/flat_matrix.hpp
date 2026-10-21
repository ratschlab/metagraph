#ifndef __FLAT_MATRIX_HPP__
#define __FLAT_MATRIX_HPP__

#include <vector>
#include <memory>

#include "annotation/binary_matrix/base/binary_matrix.hpp"

class bit_vector_sd;


namespace mtg {
namespace annot {
namespace matrix {

template <typename BitVector = bit_vector_sd>
class RowFlat : public RowMajor, public GetEntrySupport {
  public:
    RowFlat() : compressed_rows_(new BitVector()) {}
    RowFlat(const std::function<void(const RowCallback&)> &call_rows,
            uint64_t num_columns,
            uint64_t num_rows,
            uint64_t num_set_bits);

    // number of ones in the matrix
    uint64_t num_relations() const { return num_relations_; }
    uint64_t num_columns() const { return num_columns_; }
    uint64_t num_rows() const { return num_rows_; }

    bool get(Row row, Column column) const;
    SetBitPositions get_row(Row row) const;
    std::vector<Row> get_column(Column column) const;

    // Iterate the rows in order, invoking `callback` once per row with its
    // `SetBitPositions`. Mirrors `BRWT::call_rows` so RowDiff converters can
    // walk both representations through the same interface.
    void call_rows(const std::function<void(const SetBitPositions &)> &callback) const;

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    // Append a freshly-serialized RowFlat body to `filename`. The file must
    // already exist; callers write any preceding header (label encoder,
    // version magic, RowDiff anchors, …) before calling this.
    static void serialize(const std::string &filename,
                          const std::function<void(const RowCallback&)> &call_rows,
                          uint64_t num_columns,
                          uint64_t num_rows,
                          uint64_t num_set_bits);

    const BitVector& data() const { return *compressed_rows_; }

  private:
    std::unique_ptr<BitVector> compressed_rows_;
    uint64_t num_relations_ = 0;
    uint64_t num_columns_ = 0;
    uint64_t num_rows_ = 0;
};

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __FLAT_MATRIX_HPP__
