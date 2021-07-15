#ifndef __TUPLE_CSC_MATRIX_HPP__
#define __TUPLE_CSC_MATRIX_HPP__

#include <vector>

#include "annotation/int_matrix/base/int_matrix.hpp"
#include "common/logger.hpp"


namespace mtg {
namespace annot {
namespace matrix {

/**
 * Multi-Value Compressed Sparse Column Matrix (column-rank extended)
 *
 * Matrix which stores the non-empty tuples externally and indexes their
 * positions in a binary matrix. These values are indexed by rank1 called
 * on binary columns of the indexing matrix.
 */
template <class BaseMatrix,
          class Values = sdsl::int_vector<>,
          class Delims = bit_vector_smart>
class TupleCSCMatrix : public MultiIntMatrix {
  public:
    TupleCSCMatrix() {}

    TupleCSCMatrix(BaseMatrix&& index_matrix)
      : binary_matrix_(std::move(index_matrix)) {}

    // return tuple sizes (if not zero) at each entry
    RowValues get_row_values(Row row) const;

    std::vector<RowValues>
    get_row_values(const std::vector<Row> &rows) const;

    uint64_t num_attributes() const;

    // return entries of the matrix -- where each entry is a set of integers
    RowTuples get_row_tuples(Row row) const;

    std::vector<RowTuples>
    get_row_tuples(const std::vector<Row> &rows) const;

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

    bool load_tuples(std::istream &in);
    void serialize_tuples(std::ostream &out) const;

    const BaseMatrix& get_binary_matrix() const { return binary_matrix_; }

  private:
    BaseMatrix binary_matrix_;
    std::vector<Delims> delimiters_;
    std::vector<Values> column_values_;
};


template <class BaseMatrix, class Values, class Delims>
inline typename TupleCSCMatrix<BaseMatrix, Values, Delims>::RowValues
TupleCSCMatrix<BaseMatrix, Values, Delims>::get_row_values(Row row) const {
    const auto &column_ranks = binary_matrix_.get_column_ranks(row);
    RowValues row_values;
    row_values.reserve(column_ranks.size());
    for (auto [j, r] : column_ranks) {
        assert(r >= 1 && "matches can't have zero-rank");
        size_t tuple_size = delimiters_[j].select1(r + 1) - delimiters_[j].select1(r) - 1;
        row_values.emplace_back(j, tuple_size);
    }
    return row_values;
}

template <class BaseMatrix, class Values, class Delims>
inline std::vector<typename TupleCSCMatrix<BaseMatrix, Values, Delims>::RowValues>
TupleCSCMatrix<BaseMatrix, Values, Delims>::get_row_values(const std::vector<Row> &rows) const {
    const auto &column_ranks = binary_matrix_.get_column_ranks(rows);
    std::vector<RowValues> row_values(rows.size());
    // TODO: reshape?
    for (size_t i = 0; i < rows.size(); ++i) {
        row_values[i].reserve(column_ranks[i].size());
        for (auto [j, r] : column_ranks[i]) {
            assert(r >= 1 && "matches can't have zero-rank");
            size_t tuple_size = delimiters_[j].select1(r + 1) - delimiters_[j].select1(r) - 1;
            row_values[i].emplace_back(j, tuple_size);
        }
    }
    return row_values;
}

template <class BaseMatrix, class Values, class Delims>
uint64_t TupleCSCMatrix<BaseMatrix, Values, Delims>::num_attributes() const {
    uint64_t num_attributes = 0;
    for (size_t j = 0; j < column_values_.size(); ++j) {
        num_attributes += column_values_[j].size();
    }
    return num_attributes;
}

template <class BaseMatrix, class Values, class Delims>
inline typename TupleCSCMatrix<BaseMatrix, Values, Delims>::RowTuples
TupleCSCMatrix<BaseMatrix, Values, Delims>::get_row_tuples(Row row) const {
    const auto &column_ranks = binary_matrix_.get_column_ranks(row);
    RowTuples row_tuples;
    row_tuples.reserve(column_ranks.size());
    for (auto [j, r] : column_ranks) {
        assert(r >= 1 && "matches can't have zero-rank");
        size_t begin = delimiters_[j].select1(r) + 1 - r;
        size_t end = delimiters_[j].select1(r + 1) - r;
        Tuple tuple;
        tuple.reserve(end - begin);
        for (size_t t = begin; t < end; ++t) {
            tuple.push_back(column_values_[j][t]);
        }
        row_tuples.emplace_back(j, std::move(tuple));
    }
    return row_tuples;
}

template <class BaseMatrix, class Values, class Delims>
inline std::vector<typename TupleCSCMatrix<BaseMatrix, Values, Delims>::RowTuples>
TupleCSCMatrix<BaseMatrix, Values, Delims>::get_row_tuples(const std::vector<Row> &rows) const {
    const auto &column_ranks = binary_matrix_.get_column_ranks(rows);
    std::vector<RowTuples> row_tuples(rows.size());
    // TODO: reshape?
    for (size_t i = 0; i < rows.size(); ++i) {
        row_tuples[i].reserve(column_ranks[i].size());
        for (auto [j, r] : column_ranks[i]) {
            assert(r >= 1 && "matches can't have zero-rank");
            size_t begin = delimiters_[j].select1(r) + 1 - r;
            size_t end = delimiters_[j].select1(r + 1) - r;
            Tuple tuple;
            tuple.reserve(end - begin);
            for (size_t t = begin; t < end; ++t) {
                tuple.push_back(column_values_[j][t]);
            }
            row_tuples[i].emplace_back(j, std::move(tuple));
        }
    }
    return row_tuples;
}

template <class BaseMatrix, class Values, class Delims>
inline bool TupleCSCMatrix<BaseMatrix, Values, Delims>::load(std::istream &in) {
    return binary_matrix_.load(in) && load_tuples(in);
}

template <class BaseMatrix, class Values, class Delims>
inline bool TupleCSCMatrix<BaseMatrix, Values, Delims>::load_tuples(std::istream &in) {
    delimiters_.clear();
    column_values_.clear();

    delimiters_.resize(num_columns());
    column_values_.resize(num_columns());
    for (size_t j = 0; j < column_values_.size(); ++j) {
        try {
            delimiters_[j].load(in);
            column_values_[j].load(in);
        } catch (...) {
            common::logger->error("Couldn't load multi-integer values for column {}", j);
            return false;
        }
    }
    return true;
}

template <class BaseMatrix, class Values, class Delims>
inline void TupleCSCMatrix<BaseMatrix, Values, Delims>::serialize(std::ostream &out) const {
    binary_matrix_.serialize(out);
    serialize_tuples(out);
}

template <class BaseMatrix, class Values, class Delims>
inline void TupleCSCMatrix<BaseMatrix, Values, Delims>::serialize_tuples(std::ostream &out) const {
    for (size_t j = 0; j < column_values_.size(); ++j) {
        delimiters_[j].serialize(out);
        column_values_[j].serialize(out);
    }
}

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __TUPLE_CSC_MATRIX_HPP__
