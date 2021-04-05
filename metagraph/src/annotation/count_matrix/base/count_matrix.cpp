#include "count_matrix.hpp"


namespace mtg {
namespace annot {
namespace matrix {

std::vector<std::pair<CountMatrix::Column, size_t /* count */>>
CountMatrix::sum_row_counts(const std::vector<std::pair<Row, size_t>> &index_counts,
                            size_t min_count,
                            size_t count_cap) const {
    assert(count_cap >= min_count);

    if (!count_cap)
        return {};

    min_count = std::max(min_count, size_t(1));

    size_t total_sum_count = 0;
    for (const auto &pair : index_counts) {
        total_sum_count += pair.second;
    }

    if (total_sum_count < min_count)
        return {};

    std::vector<size_t> col_counts(num_columns(), 0);

    for (auto [i, count] : index_counts) {
        for (auto [j, weight] : get_row_counts(i)) {
            assert(j < col_counts.size());
            col_counts[j] += count * weight;
        }
    }

    std::vector<std::pair<uint64_t, size_t>> result;
    result.reserve(col_counts.size());

    for (size_t j = 0; j < num_columns(); ++j) {
        if (col_counts[j] >= min_count) {
            result.emplace_back(j, std::min(col_counts[j], count_cap));
        }
    }

    return result;
}

} // namespace matrix
} // namespace annot
} // namespace mtg
