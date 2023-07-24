#include "int_matrix.hpp"


namespace mtg {
namespace annot {
namespace matrix {

using Row = BinaryMatrix::Row;
using Column = BinaryMatrix::Column;


IntMatrix::RowValues
IntMatrix::sum_row_values(const std::vector<std::pair<Row, size_t>> &index_counts,
                          size_t min_count) const {
    min_count = std::max(min_count, size_t(1));

    std::vector<Row> rows;
    rows.reserve(index_counts.size());

    size_t total_sum = 0;
    for (const auto &[i, count] : index_counts) {
        total_sum += count;
        rows.push_back(i);
    }

    if (total_sum < min_count)
        return {};

    auto row_values = get_row_values(rows);

    size_t n_cols = dynamic_cast<const BinaryMatrix &>(*this).num_columns();
    Vector<std::pair<size_t, size_t>> counts(n_cols, std::make_pair(0, 0));

    for (size_t t = 0; t < index_counts.size(); ++t) {
        auto [i, count] = index_counts[t];
        for (const auto &[j, value] : row_values[t]) {
            counts[j].first += count;
            counts[j].second += count * value;
        }
    }

    RowValues result;
    result.reserve(n_cols);

    for (size_t j = 0; j < n_cols; ++j) {
        if (counts[j].first >= min_count) {
            result.emplace_back(j, counts[j].second);
        }
    }

    return result;
}

// for each row return the sizes of all non-empty tuples
std::vector<MultiIntMatrix::RowValues>
MultiIntMatrix::get_row_values(const std::vector<Row> &rows) const {
    std::vector<RowTuples> row_tuples = get_row_tuples(rows);

    std::vector<RowValues> row_values(row_tuples.size());

    for (size_t i = 0; i < row_tuples.size(); ++i) {
        row_values[i].resize(row_tuples[i].size());
        for (size_t j = 0; j < row_tuples[i].size(); ++j) {
            row_values[i][j].first = row_tuples[i][j].first;
            row_values[i][j].second = row_tuples[i][j].second.size();
        }
    }

    return row_values;
}

} // namespace matrix
} // namespace annot
} // namespace mtg
