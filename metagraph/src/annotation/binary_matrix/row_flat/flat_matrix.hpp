#ifndef __FLAT_MATRIX_HPP__
#define __FLAT_MATRIX_HPP__

#include <vector>

#include "common/vectors/bit_vector.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"


template <typename BitVector = bit_vector_sd>
class RowConcatenated : public BinaryMatrix {
  public:
    RowConcatenated() : compressed_rows_(new BitVector()) {}
    RowConcatenated(const std::function<void(const RowCallback&)> &call_rows,
                    uint64_t num_columns,
                    uint64_t num_rows,
                    uint64_t num_set_bits);

    uint64_t num_columns() const { return num_columns_; }
    uint64_t num_rows() const { return num_rows_; }

    bool get(Row row, Column column) const;
    SetBitPositions get_row(Row row) const;
    std::vector<Row> get_column(Column column) const;
    std::vector<size_t> get_column_counts() const;
    size_t get_column_count(Column column) const;

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    // number of ones in the matrix
    uint64_t num_relations() const { return compressed_rows_->num_set_bits(); }

    const BitVector& data() const { return *compressed_rows_; }

  private:
    std::unique_ptr<BitVector> compressed_rows_;
    uint64_t num_columns_ = 0;
    uint64_t num_rows_ = 0;
};

#endif // __FLAT_MATRIX_HPP__
