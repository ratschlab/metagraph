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

IntMatrix::SetBitPositions IntMatrix::get_row(Row i) const {
    RowValues row = get_row_values(i);
    SetBitPositions result(row.size());
    for (size_t k = 0; k < row.size(); ++k) {
        result[k] = row[k].first;
    }
    return result;
}

std::vector<IntMatrix::SetBitPositions>
IntMatrix::get_rows(const std::vector<Row> &row_ids) const {
    std::vector<SetBitPositions> result;
    result.reserve(row_ids.size());

    for (auto&& row : get_row_values(row_ids)) {
        result.emplace_back(row.size());
        for (size_t k = 0; k < row.size(); ++k) {
            result.back()[k] = row[k].first;
        }
        row = RowValues();
    }

    return result;
}

IntMatrix::Value IntMatrix::get_value(Row row, Column column) const {
    RowValues values = get_row_values(row);

    for (size_t k = 0; k < values.size(); ++k) {
        if (values[k].first == column)
            return values[k].second;
    }

    return 0;
}

// return the positions of all non-empty tuples in the row
MultiIntMatrix::SetBitPositions MultiIntMatrix::get_row(Row i) const {
    RowTuples row = get_row_tuples(i);
    SetBitPositions result(row.size());
    for (size_t k = 0; k < row.size(); ++k) {
        result[k] = row[k].first;
    }
    return result;
}

// for each row return the positions of all non-empty tuples
std::vector<MultiIntMatrix::SetBitPositions>
MultiIntMatrix::get_rows(const std::vector<Row> &row_ids) const {
    std::vector<SetBitPositions> result;
    result.reserve(row_ids.size());

    for (auto&& row : get_row_tuples(row_ids)) {
        result.emplace_back(row.size());
        for (size_t k = 0; k < row.size(); ++k) {
            result.back()[k] = row[k].first;
        }
        row = RowTuples();
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

MultiIntMatrix::Tuple MultiIntMatrix::get_tuple(Row row, Column column) const {
    RowTuples row_tuples = get_row_tuples(row);
    for (size_t i = 0; i < row_tuples.size(); ++i) {
        if (row_tuples[i].first == column)
            return row_tuples[i].second;
    }

    return {};
}

} // namespace matrix
} // namespace annot
} // namespace mtg
