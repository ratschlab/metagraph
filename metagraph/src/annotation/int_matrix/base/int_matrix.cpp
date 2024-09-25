#include "int_matrix.hpp"

#include <progress_bar.hpp>

#include "common/logger.hpp"
#include "common/threads/threading.hpp"
#include "common/vectors/vector_algorithm.hpp"


namespace mtg {
namespace annot {
namespace matrix {

using Row = BinaryMatrix::Row;
using Column = BinaryMatrix::Column;

void IntMatrix::call_row_values(const std::function<void(uint64_t, const RowValues&, size_t)> &callback,
                                bool ordered) const {
    size_t n = get_binary_matrix().num_rows();
    size_t batch_size = (n + get_num_threads() - 1) / get_num_threads();
    size_t rows_per_update = 10000;
    ProgressBar progress_bar(n, "Streaming rows", std::cerr, !common::get_verbose());

    std::ignore = ordered;

    #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
    for (size_t k = 0; k < n; k += batch_size) {
        size_t begin = k;
        size_t end = std::min(begin + batch_size, n);
        size_t thread_id = k / batch_size;
        std::vector<Row> row_ids;
        for (size_t row_batch_i = begin; row_batch_i < end; row_batch_i += rows_per_update) {
            size_t row_batch_end = std::min(row_batch_i + rows_per_update, end);
            row_ids.resize(row_batch_end - row_batch_i);
            std::iota(row_ids.begin(), row_ids.end(), row_batch_i);
            auto row_vals = get_row_values(row_ids);
            for (size_t i = 0; i < row_ids.size(); ++i) {
                callback(i + row_batch_i, row_vals[i], thread_id);
            }
            progress_bar += row_ids.size();
        }
    }
}

std::vector<VectorMap<uint64_t, size_t>> IntMatrix::get_histograms(const std::vector<size_t> &min_counts,
                                                                   sdsl::bit_vector *unmark_discarded) const {
    common::logger->trace("Calculating count histograms");
    bool parallel = get_num_threads() > 0;
    std::vector<VectorMap<uint64_t, size_t>> hists_map(get_binary_matrix().num_columns());
    std::atomic_thread_fence(std::memory_order_release);
    call_row_values([&](uint64_t row_i, const auto &row, size_t) {
        if (min_counts.size()) {
            bool keep_row = false;
            for (const auto &[j, c] : row) {
                if (c >= min_counts[j]) {
                    keep_row = true;
                    break;
                }
            }

            if (!keep_row) {
                if (unmark_discarded)
                    unset_bit(unmark_discarded->data(), row_i, parallel, std::memory_order_relaxed);

                return;
            }
        }

        Vector<uint64_t> counts(hists_map.size());
        for (const auto &[j, raw_c] : row) {
            counts[j] = raw_c;
        }

        #pragma omp critical
        {
            for (size_t j = 0; j < counts.size(); ++j) {
                ++hists_map[j][counts[j]];
            }
        }
    }, false);
    std::atomic_thread_fence(std::memory_order_acquire);

    return hists_map;
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
