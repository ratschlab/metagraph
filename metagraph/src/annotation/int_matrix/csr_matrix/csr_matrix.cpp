#include "csr_matrix.hpp"

#include <algorithm>
#include <numeric>

#include "common/algorithms.hpp"
#include "common/utils/template_utils.hpp"


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
CSRMatrix::get_row_values(const std::vector<Row> &rows, size_t num_threads) const {
    std::vector<RowValues> row_values(rows.size());
    #pragma omp parallel for num_threads(num_threads)
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

CSRMatrix::SetBitPositions CSRMatrix::get_row(Row row) const {
    assert(row < vector_.size());
    return utils::get_firsts<SetBitPositions>(vector_[row]);
}

std::vector<std::pair<CSRMatrix::Column, size_t /* count */>>
CSRMatrix::sum_rows(const std::vector<std::pair<Row, size_t>> &index_counts, size_t min_count) const {
    if (index_counts.empty())
        return {};

    min_count = std::max<size_t>(min_count, 1);

    size_t total_sum_count = 0;
    for (const auto &[row, count] : index_counts) {
        total_sum_count += count;
    }
    if (total_sum_count < min_count)
        return {};

    auto call_bits = [&](const auto &callback) {
        for (const auto &[row, count] : index_counts) {
            for (const auto &j : vector_[row]) {
                callback(utils::get_first(j), count);
            }
        }
    };
    return utils::accumulate_counts(call_bits, num_columns(), min_count);
}

std::vector<CSRMatrix::Row> CSRMatrix::get_column(Column column) const {
    std::vector<Row> result;
    for (uint64_t i = 0; i < vector_.size(); ++i) {
        const auto &row = vector_[i];
        if (std::find_if(row.begin(), row.end(), [&](const auto &p) { return p.first == column; }) != row.end()) {
            result.push_back(i);
        }
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
