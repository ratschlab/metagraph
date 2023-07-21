#ifndef __RAINBOWFISH_HPP__
#define __RAINBOWFISH_HPP__

#include <vector>

#include "common/vectors/bit_vector_sdsl.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "annotation/binary_matrix/row_flat/flat_matrix.hpp"


namespace mtg {
namespace annot {
namespace matrix {

// TODO: remove (can be replaced with Rainbow<RowFlat>)
class Rainbowfish : public RainbowMatrix {
  public:
    Rainbowfish() {}
    Rainbowfish(const std::function<void(RowCallback)> &call_rows,
                uint64_t num_columns,
                uint64_t buffer_size = static_cast<uint64_t>(-1));

    uint64_t num_columns() const { return num_columns_; }
    uint64_t num_rows() const;
    uint64_t num_distinct_rows() const;

    // row is in [0, num_rows), column is in [0, num_columns)
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
    using ReducedMatrixType = RowFlat<bit_vector_small>;

    uint64_t num_columns_ = 0;
    uint64_t num_relations_ = 0;
    uint64_t buffer_size_ = static_cast<uint64_t>(-1);

    bit_vector_rrr<> row_codes_;
    bit_vector_rrr<> row_code_delimiters_;
    std::vector<std::unique_ptr<BinaryMatrix>> reduced_matrix_;

    uint64_t get_code(Row row) const;
    std::vector<SetBitPositions> codes_to_rows(const std::vector<uint64_t> &rows) const {
        std::vector<SetBitPositions> result;
        result.reserve(rows.size());
        for (uint64_t c : rows) {
            result.push_back(reduced_matrix_[c / buffer_size_]->get_rows({ c % buffer_size_ })[0]);
        }
        return result;
    }
};

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __RAINBOWFISH_HPP__
