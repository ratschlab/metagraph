#ifndef __RAINBOW_HPP__
#define __RAINBOW_HPP__

#include <vector>

#include "common/vectors/bit_vector_sdsl.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"


namespace mtg {
namespace annot {
namespace binmat {

template <class MatrixType>
class Rainbow : public RainbowMatrix {
  public:
    typedef MatrixType matrix_type;

    Rainbow() {}

    Rainbow(MatrixType&& reduced_matrix,
            sdsl::bit_vector&& row_codes_,
            bit_vector_rrr<>&& row_code_delimiters_,
            uint64_t num_relations);

    // using CallColumn = std::function<void(const std::unique_ptr<bit_vector> &)>;
    // Rainbow(const std::function<void(const CallColumn &)> &get_columns);

    uint64_t num_columns() const { return reduced_matrix_.num_columns(); }
    uint64_t num_rows() const;
    uint64_t num_distinct_rows() const { return reduced_matrix_.num_rows(); }

    // row is in [0, num_rows), column is in [0, num_columns)
    bool get(Row row, Column column) const;
    SetBitPositions get_row(Row row) const;
    std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const;
    // Return unique rows (in arbitrary order) and update the row indexes
    // in |rows| to point to their respective rows in the vector returned.
    std::vector<SetBitPositions> get_rows(std::vector<Row> *rows,
                                          size_t num_threads = 1) const;
    std::vector<Row> get_column(Column column) const;

    void slice_columns(const std::vector<Column> &columns,
                       const ColumnCallback &callback) const;

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    // number of ones in the matrix
    uint64_t num_relations() const { return num_relations_; }

    const MatrixType& get_reduced_matrix() const { return reduced_matrix_; }

  private:
    uint64_t num_relations_ = 0;

    sdsl::bit_vector row_codes_;
    bit_vector_rrr<> row_code_delimiters_;
    MatrixType reduced_matrix_;

    uint64_t get_code(Row row) const;
};

} // namespace binmat
} // namespace annot
} // namespace mtg

#endif // __RAINBOW_HPP__
