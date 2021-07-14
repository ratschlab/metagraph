#include "int_matrix.hpp"


namespace mtg {
namespace annot {
namespace matrix {

IntMatrix::RowValues
IntMatrix::sum_row_values(const std::vector<std::pair<Row, size_t>> &index_counts,
                          size_t min,
                          size_t cap) const {
    assert(cap >= min);

    if (!cap)
        return {};

    min = std::max(min, size_t(1));

    size_t total_sum = 0;
    for (const auto &pair : index_counts) {
        total_sum += pair.second;
    }

    if (total_sum < min)
        return {};

    std::vector<size_t> sum_row(num_columns(), 0);

    for (auto [i, count] : index_counts) {
        for (auto [j, value] : get_row_values(i)) {
            assert(j < sum_row.size());
            sum_row[j] += count * value;
        }
    }

    RowValues result;
    result.reserve(sum_row.size());

    for (size_t j = 0; j < num_columns(); ++j) {
        if (sum_row[j] >= min) {
            result.emplace_back(j, std::min(sum_row[j], cap));
        }
    }

    return result;
}


// return tuple sizes (if not zero) at each entry
MultiIntMatrix::RowValues MultiIntMatrix::get_row_values(Row row) const {
    RowTuples row_tuples = get_row_tuples(row);

    RowValues row_values(row_tuples.size());

    for (size_t i = 0; i < row_tuples.size(); ++i) {
        row_values[i].first = row_tuples[i].first;
        row_values[i].second = row_tuples[i].second.size();
    }

    return row_values;
}

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
