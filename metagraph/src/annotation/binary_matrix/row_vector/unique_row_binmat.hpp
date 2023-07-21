#ifndef __UNIQUE_ROW_BINMAT_HPP__
#define __UNIQUE_ROW_BINMAT_HPP__

#include <vector>

#include "annotation/binary_matrix/base/binary_matrix.hpp"


namespace mtg {
namespace annot {
namespace matrix {

class UniqueRowBinmat : public RainbowMatrix {
  public:
    explicit UniqueRowBinmat(uint64_t num_rows = 0);

    UniqueRowBinmat(std::vector<SetBitPositions>&& unique_rows,
                    std::vector<uint32_t>&& row_rank,
                    uint32_t num_columns);

    UniqueRowBinmat(const std::function<void(const RowCallback &)> &call_rows,
                    uint32_t num_columns);

    uint64_t num_columns() const { return num_columns_; }
    uint64_t num_rows() const { return row_rank_.size(); }
    uint64_t num_distinct_rows() const { return unique_rows_.size(); }

    std::vector<Row> get_column(Column column) const;

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    // number of ones in the matrix
    uint64_t num_relations() const { return num_relations_; }
    // matrix density
    double density() const;

  private:
    uint64_t get_code(Row row) const { return row_rank_[row]; }
    std::vector<SetBitPositions> codes_to_rows(const std::vector<uint64_t> &rows) const {
        std::vector<SetBitPositions> result(rows.size());
        for (size_t i = 0; i < rows.size(); ++i) {
            result[i] = unique_rows_[rows[i]];
        }
        return result;
    }

    uint32_t num_columns_ = 0;
    uint32_t num_relations_ = 0;
    std::vector<SetBitPositions> unique_rows_;
    std::vector<uint32_t> row_rank_;
};

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __UNIQUE_ROW_BINMAT_HPP__
