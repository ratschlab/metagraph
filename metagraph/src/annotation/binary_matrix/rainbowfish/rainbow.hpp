#ifndef __RAINBOW_HPP__
#define __RAINBOW_HPP__

#include <vector>

#include "common/vectors/bit_vector_sdsl.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"


namespace mtg {
namespace annot {
namespace matrix {

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

    uint64_t num_columns() const override { return reduced_matrix_.num_columns(); }
    uint64_t num_rows() const override;
    uint64_t num_distinct_rows() const override { return reduced_matrix_.num_rows(); }

    std::vector<Row> get_column(Column column) const override;

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    // number of ones in the matrix
    uint64_t num_relations() const override { return num_relations_; }

    const MatrixType& get_reduced_matrix() const { return reduced_matrix_; }

  private:
    uint64_t num_relations_ = 0;

    // TODO: try out sdsl::dac_vector<> or other compression methods
    sdsl::bit_vector row_codes_;
    bit_vector_rrr<> row_code_delimiters_;
    MatrixType reduced_matrix_;

    uint64_t get_code(Row row) const override;
    std::vector<bit_vector_smart> codes_to_rows(const std::vector<uint64_t> &rows) const override {
        auto ext_rows = reduced_matrix_.get_rows(rows);
        std::vector<bit_vector_smart> result;
        result.reserve(ext_rows.size());
        for (const auto &row : ext_rows) {
            result.emplace_back(
                [&](const auto &callback) {
                    for (auto j : row) {
                        callback(j);
                    }
                },
                num_columns(),
                row.size()
            );
        }

        return result;
    }
};

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __RAINBOW_HPP__
