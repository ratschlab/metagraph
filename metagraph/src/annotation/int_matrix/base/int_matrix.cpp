#include "int_matrix.hpp"


namespace mtg {
namespace annot {
namespace matrix {

std::vector<IntMatrix::RowValueDiffs>
IntMatrix::get_row_value_diffs(const std::vector<Row> &rows) const {
    if (rows.empty())
        return {};

    auto row_values = get_row_values(rows);
    std::vector<RowValueDiffs> row_value_diffs(rows.size());
    std::sort(row_values[0].begin(), row_values[0].end());

    for (size_t i = 1; i < rows.size(); ++i) {
        auto &result = row_value_diffs[i];
        std::sort(row_values[i].begin(), row_values[i].end());
        auto it = row_values[i].begin();
        auto it2 = row_values[i - 1].begin();
        while (it != row_values[i].end() && it2 != row_values[i - 1].end()) {
            if (it->first < it2->first) {
                result.push_back(*it);
                ++it;
            } else if (it->first > it2->first) {
                result.emplace_back(it2->first, -static_cast<ValueDiff>(it2->second));
                ++it2;
            } else {
                if (int64_t diff = it->second - it2->second)
                    result.emplace_back(it->first, diff);
                ++it;
                ++it2;
            }
        }
        std::copy(it, row_values[i].end(), std::back_inserter(result));
        std::transform(it2, row_values[i - 1].end(), std::back_inserter(result),
                       [](const auto &val) {
                           return std::make_pair(val.first,
                                                 -static_cast<ValueDiff>(val.second));
                       });
    }

    row_value_diffs[0] = RowValueDiffs(row_values[0].begin(), row_values[0].end());

    return row_value_diffs;
}

MultiIntMatrix::RowTuples MultiIntMatrix::get_row_tuple_diff(Row a, Row b) const {
    auto row_tuples = get_row_tuples(std::vector<Row>{ a, b });
    for (auto &[j, tuple] : row_tuples[0]) {
        for (auto &c : tuple) {
            c += 1;
        }
    }
    std::sort(row_tuples[0].begin(), row_tuples[0].end());
    std::sort(row_tuples[1].begin(), row_tuples[1].end());

    RowTuples result;
    auto it = row_tuples[1].begin();
    auto it2 = row_tuples[0].begin();
    while (it != row_tuples[1].end() && it2 != row_tuples[0].end()) {
        if (it->first < it2->first) {
            assert(std::is_sorted(it->second.begin(), it->second.end()));
            result.push_back(*it);
            ++it;
        } else if (it->first > it2->first) {
            assert(std::is_sorted(it2->second.begin(), it2->second.end()));
            result.push_back(*it2);
            ++it2;
        } else {
            assert(std::is_sorted(it->second.begin(), it->second.end()));
            assert(std::is_sorted(it2->second.begin(), it2->second.end()));
            result.emplace_back(it->first, Tuple{});
            std::set_symmetric_difference(it->second.begin(), it->second.end(),
                                          it2->second.begin(), it2->second.end(),
                                          std::back_inserter(result.back().second));
            if (result.back().second.empty())
                result.pop_back();

            ++it;
            ++it2;
        }
    }
    std::copy(it, row_tuples[1].end(), std::back_inserter(result));
    std::copy(it2, row_tuples[0].end(), std::back_inserter(result));
    return result;
}

std::vector<MultiIntMatrix::RowTuples>
MultiIntMatrix::get_row_tuple_diffs(const std::vector<Row> &rows, const RowTuples *first_tuple) const {
    if (rows.empty())
        return {};

    if (rows.size() == 1) {
        if (first_tuple)
            return { *first_tuple };

        return get_row_tuples(rows);
    }

    auto row_tuples = get_row_tuples(rows);
    std::vector<RowTuples> row_tuple_diffs(rows.size() - 1);
    std::sort(row_tuples[0].begin(), row_tuples[0].end());
    for (size_t i = 1; i < rows.size(); ++i) {
        auto &result = row_tuple_diffs[i];
        std::sort(row_tuples[i].begin(), row_tuples[i].end());
        auto it = row_tuples[i].begin();
        auto it2 = row_tuples[i - 1].begin();
        while (it != row_tuples[i].end() && it2 != row_tuples[i - 1].end()) {
            if (it->first < it2->first) {
                result.push_back(*it);
                ++it;
            } else if (it->first > it2->first) {
                result.push_back(*it2);
                ++it2;
            } else {
                result.emplace_back(it->first, Tuple{});
                std::set_symmetric_difference(it->second.begin(), it->second.end(),
                                              it2->second.begin(), it2->second.end(),
                                              std::back_inserter(result.back().second));
                ++it;
                ++it2;
            }
        }
        std::copy(it, row_tuples[i].end(), std::back_inserter(result));
        std::copy(it2, row_tuples[i - 1].end(), std::back_inserter(result));
    }
    row_tuple_diffs[0] = std::move(row_tuples[0]);

    return row_tuple_diffs;
}

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
