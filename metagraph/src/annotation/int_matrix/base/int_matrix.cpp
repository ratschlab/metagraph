#include "int_matrix.hpp"

#include <tsl/hopscotch_map.h>

#include "common/algorithms.hpp"


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

    auto sum_rows = [&](auto &counts) {
        for (size_t t = 0; t < index_counts.size(); ++t) {
            auto [i, count] = index_counts[t];
            for (const auto &[j, value] : row_values[t]) {
                counts[j].first += count;
                counts[j].second += count * value;
            }
        }
    };

    size_t num_cols = dynamic_cast<const BinaryMatrix &>(*this).num_columns();
    constexpr size_t kMinReserve = 1024;
    if (num_cols < utils::kDenseCountThreshold) {
        // For few columns, counting with a dense vector is faster than a hash table
        Vector<std::pair<size_t, size_t>> counts(num_cols);
        sum_rows(counts);

        RowValues result;
        result.reserve(std::min<size_t>(num_cols, kMinReserve));

        for (size_t j = 0; j < num_cols; ++j) {
            if (counts[j].first >= min_count) {
                result.emplace_back(j, counts[j].second);
            }
        }
        return result;
    }

    tsl::hopscotch_map<Column, std::pair<size_t, size_t>> counts;
    counts.reserve(std::min<size_t>(num_cols, kMinReserve));
    sum_rows(counts);

    RowValues result;
    result.reserve(std::min<size_t>(num_cols, kMinReserve));

    for (const auto &[j, pair] : counts) {
        if (pair.first >= min_count) {
            result.emplace_back(j, pair.second);
        }
    }
    return result;
}

} // namespace matrix
} // namespace annot
} // namespace mtg
