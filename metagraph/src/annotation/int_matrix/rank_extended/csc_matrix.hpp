#ifndef __CSC_MATRIX_HPP__
#define __CSC_MATRIX_HPP__

#include <vector>

#include "annotation/int_matrix/base/int_matrix.hpp"
#include "common/logger.hpp"


namespace mtg {
namespace annot {
namespace matrix {

/**
 * Compressed Sparse Column Matrix (column-rank extended)
 *
 * Matrix which stores the non-zero values externally and indexes their
 * positions in a binary matrix. These values are indexed by rank1 called
 * on binary columns of the indexing matrix.
 */
template <class BaseMatrix, class ColumnValues = sdsl::int_vector<>>
class CSCMatrix : public BinaryMatrix, public IntMatrix, public GetEntrySupport {
  public:
    CSCMatrix() {}

    CSCMatrix(BaseMatrix&& index_matrix,
              std::vector<ColumnValues>&& column_values)
      : binary_matrix_(std::move(index_matrix)),
        column_values_(column_values) {}

    // row is in [0, num_rows), column is in [0, num_columns)
    std::vector<RowValues> get_row_values(const std::vector<Row> &rows) const;

    uint64_t num_columns() const { return binary_matrix_.num_columns(); }
    uint64_t num_rows() const { return binary_matrix_.num_rows(); }
    uint64_t num_relations() const { return binary_matrix_.num_relations(); }

    // row is in [0, num_rows), column is in [0, num_columns)
    bool get(Row row, Column column) const { return binary_matrix_.get(row, column); }
    std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const {
        return binary_matrix_.get_rows(rows);
    }
    std::vector<Row> get_column(Column column) const {
        return binary_matrix_.get_column(column);
    }

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    const BaseMatrix& get_binary_matrix() const { return binary_matrix_; }

  private:
    BaseMatrix binary_matrix_;
    std::vector<ColumnValues> column_values_;
};


template <class BaseMatrix, class ColumnValues>
inline std::vector<typename CSCMatrix<BaseMatrix, ColumnValues>::RowValues>
CSCMatrix<BaseMatrix, ColumnValues>::get_row_values(const std::vector<Row> &rows) const {
    const auto &column_ranks = binary_matrix_.get_column_ranks(rows);
    std::vector<RowValues> row_values(rows.size());
    // TODO: reshape?
    for (size_t i = 0; i < rows.size(); ++i) {
        row_values[i].reserve(column_ranks[i].size());
        for (auto [j, r] : column_ranks[i]) {
            assert(r >= 1 && "matches can't have zero-rank");
            row_values[i].emplace_back(j, column_values_[j][r - 1]);
        }
    }
    return row_values;
}

template <class BaseMatrix, class ColumnValues>
inline bool CSCMatrix<BaseMatrix, ColumnValues>::load(std::istream &in) {
    column_values_.clear();

    if (!binary_matrix_.load(in))
        return false;

    column_values_.resize(num_columns());
    for (size_t j = 0; j < column_values_.size(); ++j) {
        try {
            column_values_[j].load(in);
        } catch (...) {
            common::logger->error("Couldn't load integer values for column {}", j);
            return false;
        }
    }
    return true;
}

template <class BaseMatrix, class ColumnValues>
inline void CSCMatrix<BaseMatrix, ColumnValues>::serialize(std::ostream &out) const {
    binary_matrix_.serialize(out);
    for (const auto &col : column_values_) {
        col.serialize(out);
    }
}

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __CSC_MATRIX_HPP__
