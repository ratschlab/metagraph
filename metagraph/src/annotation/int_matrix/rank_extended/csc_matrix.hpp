#ifndef __CSC_MATRIX_HPP__
#define __CSC_MATRIX_HPP__

#include <algorithm>
#include <utility>
#include <vector>

#include "annotation/binary_matrix/base/decode_budget.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "common/logger.hpp"
#include "common/utils/template_utils.hpp"


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
        column_values_(std::move(column_values)) {}

    // row is in [0, num_rows), column is in [0, num_columns)
    std::vector<RowValues> get_row_values(std::span<const Row> rows,
                                          size_t num_threads = 1) const;

    // The budget-aware decode path (decode_budget.hpp, as TupleCSCMatrix::row_tuples): the
    // values of ONE row into the empty |*out|, sorted by column, every byte charged to
    // |budget| before it is allocated — the column ranks through |scratch| (the index
    // matrix's row_column_ranks()), then the row's vector at its exact size. |*bytes|
    // receives what is charged for |*out| (released by the caller when it frees the row).
    // false: a charge did not fit. Only for index matrices with row_column_ranks() (BRWT,
    // ColumnMajor): a member of a class template is instantiated only where it is used.
    bool row_values(Row row, ChargedBuffer<std::pair<Column, uint64_t>> *scratch,
                    RowValues *out, uint64_t *bytes, DecodeBudget &budget) const;

    uint64_t num_columns() const { return binary_matrix_.num_columns(); }
    uint64_t num_rows() const { return binary_matrix_.num_rows(); }
    uint64_t num_relations() const { return binary_matrix_.num_relations(); }

    // row is in [0, num_rows), column is in [0, num_columns)
    bool get(Row row, Column column) const { return binary_matrix_.get(row, column); }
    std::vector<SetBitPositions> get_rows(std::span<const Row> rows) const {
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
CSCMatrix<BaseMatrix, ColumnValues>::get_row_values(std::span<const Row> rows,
                                                    size_t num_threads) const {
    const auto &column_ranks = binary_matrix_.get_column_ranks(rows, num_threads);
    std::vector<RowValues> row_values(rows.size());

    #pragma omp parallel for num_threads(num_threads)
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
inline bool CSCMatrix<BaseMatrix, ColumnValues>
::row_values(Row row, ChargedBuffer<std::pair<Column, uint64_t>> *scratch,
             RowValues *out, uint64_t *bytes, DecodeBudget &budget) const {
    assert(out->empty());
    *bytes = 0;
    scratch->data.clear();
    if (!binary_matrix_.row_column_ranks(row, scratch, budget))
        return false;
    // the row is used sorted by column (as call_rows() sorts what get_row_values returns)
    std::sort(scratch->data.begin(), scratch->data.end(), utils::LessFirst());
    const uint64_t b = buffer_bytes(scratch->data.size(), sizeof((*out)[0]));
    if (!budget.charge(b))
        return false;
    *bytes = b;
    out->reserve(scratch->data.size());
    for (auto [j, r] : scratch->data) {
        assert(r >= 1 && "matches can't have zero-rank");
        out->emplace_back(j, column_values_[j][r - 1]);
    }
    return true;
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
