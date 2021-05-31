#include "csr_matrix.hpp"

#include "common/logger.hpp"


namespace mtg {
namespace annot {
namespace matrix {

CSRMatrix::CSRMatrix(Vector<RowValues>&& rows, uint64_t num_columns)
      : num_columns_(num_columns), vector_(std::move(rows)) {
    // make sure there are no columns with indexes greater than num_labels
    assert(std::all_of(vector_.begin(), vector_.end(), [&](const auto &row) {
        return std::all_of(row.begin(), row.end(),
                           [num_columns](auto pair) { return pair.first < num_columns; });
    }));
}

std::vector<CSRMatrix::RowValues>
CSRMatrix::get_row_values(const std::vector<Row> &rows) const {
    std::vector<RowValues> row_values(rows.size());
    for (size_t i = 0; i < rows.size(); ++i) {
        row_values[i] = vector_[rows[i]];
    }
    return row_values;
}

// number of non-zero entries in the matrix
uint64_t CSRMatrix::num_relations() const {
    return std::accumulate(
        vector_.begin(), vector_.end(), (uint64_t)0,
        [](uint64_t sum, const auto &v) { return sum + v.size(); }
    );
}

bool CSRMatrix::get(Row row, Column column) const {
    assert(row < vector_.size());
    return std::find_if(vector_[row].begin(), vector_[row].end(),
                        [&](const auto &pair) { return pair.first == column; })
                != vector_[row].end();
}

CSRMatrix::SetBitPositions CSRMatrix::get_row(Row row) const {
    assert(row < vector_.size());
    SetBitPositions result;
    result.reserve(vector_[row].size());
    for (const auto &[j, _] : vector_[row]) {
        result.push_back(j);
    }
    return result;
}

std::vector<CSRMatrix::Row> CSRMatrix::get_column(Column column) const {
    std::vector<Row> result;
    for (uint64_t i = 0; i < vector_.size(); ++i) {
        if (get(i, column))
            result.push_back(i);
    }
    return result;
}

bool CSRMatrix::load(std::istream &) {
    throw std::runtime_error("Not implemented");
}

void CSRMatrix::serialize(std::ostream &) const {
    throw std::runtime_error("Not implemented");
}

} // namespace matrix
} // namespace annot
} // namespace mtg
