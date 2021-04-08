#ifndef __COLUMN_RANK_MATRIX_HPP__
#define __COLUMN_RANK_MATRIX_HPP__

#include <vector>

#include "annotation/count_matrix/base/count_matrix.hpp"
#include "common/logger.hpp"


namespace mtg {
namespace annot {
namespace matrix {

/**
 * Column Rank-Extended Matrix
 *
 * Matrix which stores the non-zero values externally and indexes their
 * positions in a binary matrix. These values are indexed by rank1 called
 * on binary columns of the indexing matrix.
 */
template <class BaseMatrix>
class ColumnRankMatrix : public CountMatrix {
  public:
    ColumnRankMatrix() {}

    ColumnRankMatrix(BaseMatrix&& index_matrix,
                     std::vector<sdsl::int_vector<>>&& count_columns)
      : binary_matrix_(std::move(index_matrix)),
        count_columns_(count_columns) {}

    // row is in [0, num_rows), column is in [0, num_columns)
    ColumnCounts get_row_counts(Row row) const;

    std::vector<ColumnCounts>
    get_row_counts(const std::vector<Row> &rows) const;

    uint64_t num_columns() const { return binary_matrix_.num_columns(); }
    uint64_t num_rows() const { return binary_matrix_.num_rows(); }
    uint64_t num_relations() const { return binary_matrix_.num_relations(); }

    // row is in [0, num_rows), column is in [0, num_columns)
    bool get(Row row, Column column) const { return binary_matrix_.get(row, column); }
    SetBitPositions get_row(Row row) const { return binary_matrix_.get_row(row); }
    std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const {
        return binary_matrix_.get_rows(rows);
    }
    std::vector<Row> get_column(Column column) const {
        return binary_matrix_.get_column(column);
    }
    // get all selected rows appended with -1 and concatenated
    std::vector<Column> slice_rows(const std::vector<Row> &rows) const {
        return binary_matrix_.slice_rows(rows);
    }

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

  private:
    BaseMatrix binary_matrix_;
    std::vector<sdsl::int_vector<>> count_columns_;
};


template <class BaseMatrix>
inline typename ColumnRankMatrix<BaseMatrix>::ColumnCounts
ColumnRankMatrix<BaseMatrix>::get_row_counts(Row row) const {
    const auto &column_ranks = binary_matrix_.get_column_ranks(row);
    ColumnCounts row_counts;
    row_counts.reserve(column_ranks.size());
    for (auto [j, r] : column_ranks) {
        assert(r >= 1 && "matches can't have zero-rank");
        row_counts.emplace_back(j, count_columns_[j][r - 1]);
    }
    return row_counts;
}

template <class BaseMatrix>
inline std::vector<typename ColumnRankMatrix<BaseMatrix>::ColumnCounts>
ColumnRankMatrix<BaseMatrix>::get_row_counts(const std::vector<Row> &rows) const {
    const auto &column_ranks = binary_matrix_.get_column_ranks(rows);
    std::vector<ColumnCounts> row_counts(rows.size());
    // TODO: reshape?
    for (size_t i = 0; i < rows.size(); ++i) {
        row_counts[i].reserve(column_ranks[i].size());
        for (auto [j, r] : column_ranks[i]) {
            assert(r >= 1 && "matches can't have zero-rank");
            row_counts[i].emplace_back(j, count_columns_[j][r - 1]);
        }
    }
    return row_counts;
}

template <class BaseMatrix>
inline bool ColumnRankMatrix<BaseMatrix>::load(std::istream &in) {
    count_columns_.clear();

    if (!binary_matrix_.load(in))
        return false;

    count_columns_.resize(num_columns());
    for (size_t j = 0; j < count_columns_.size(); ++j) {
        try {
            count_columns_[j].load(in);
        } catch (...) {
            common::logger->error("Couldn't load counts for column {}", j);
            return false;
        }
    }
    return true;
}

template <class BaseMatrix>
inline void ColumnRankMatrix<BaseMatrix>::serialize(std::ostream &out) const {
    binary_matrix_.serialize(out);
    for (const auto &col : count_columns_) {
        col.serialize(out);
    }
}

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __COLUMN_RANK_MATRIX_HPP__
