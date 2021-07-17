#include "int_matrix.hpp"


namespace mtg {
namespace annot {
namespace matrix {

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

    std::vector<size_t> sum_row(num_columns(), 0);
    std::vector<size_t> counts(num_columns(), 0);

    auto row_values = get_row_values(rows);

    for (size_t t = 0; t < index_counts.size(); ++t) {
        auto [i, count] = index_counts[t];
        for (const auto &[j, value] : row_values[t]) {
            assert(j < sum_row.size());
            sum_row[j] += count * value;
            counts[j] += count;
        }
    }

    RowValues result;
    result.reserve(sum_row.size());

    for (size_t j = 0; j < num_columns(); ++j) {
        if (counts[j] >= min_count) {
            result.emplace_back(j, sum_row[j]);
        }
    }

    return result;
}


// return sizes of all non-empty tuples in the row
MultiIntMatrix::RowValues MultiIntMatrix::get_row_values(Row row) const {
    RowTuples row_tuples = get_row_tuples(row);

    RowValues row_values(row_tuples.size());

    for (size_t i = 0; i < row_tuples.size(); ++i) {
        row_values[i].first = row_tuples[i].first;
        row_values[i].second = row_tuples[i].second.size();
    }

    return row_values;
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
