#ifndef __BIT_VECTOR_UTILS_HPP__
#define __BIT_VECTOR_UTILS_HPP__

#include <vector>
#include <memory>

#include <progress_bar.hpp>
#include <sdsl/int_vector_buffer.hpp>

#include "bit_vector_sdsl.hpp"
#include "common/logger.hpp"


namespace utils {

const size_t kNumRowsInBlock = 50'000;

// Transpose all columns and call all rows.
// If there are no columns, call `num_rows` empty rows.
template <class BitmapPtr, typename Vector = std::vector<uint64_t>>
void call_rows(const std::vector<BitmapPtr> &columns,
               const std::function<void(const Vector &)> &callback,
               uint64_t num_rows = 0,
               bool show_progress = mtg::common::get_verbose()) {
    if (!columns.size()) {
        ProgressBar progress_bar(num_rows, "Constructed rows", std::cerr, !show_progress);
        Vector empty;
        while (num_rows--) {
            callback(empty);
        }
        return;
    }

    num_rows = columns[0]->size();
    std::vector<Vector> rows;

    ProgressBar progress_bar(num_rows, "Constructed rows", std::cerr, !show_progress);

    #pragma omp parallel for ordered num_threads(get_num_threads()) schedule(dynamic) private(rows)
    for (uint64_t i = 0; i < num_rows; i += kNumRowsInBlock) {
        uint64_t begin = i;
        uint64_t end = std::min(i + kNumRowsInBlock, num_rows);
        assert(begin <= end);

        rows.resize(end - begin);

        for (size_t j = 0; j < columns.size(); ++j) {
            columns[j]->call_ones_in_range(begin, end,
                [&](uint64_t i) { rows[i - begin].push_back(j); }
            );
        }
        #pragma omp ordered
        {
            for (auto &row : rows) {
                callback(row);
                ++progress_bar;
                row.resize(0);
            }
        }
    }
}

// Transpose all columns and call all rows.
// If there are no columns, call `num_rows` empty rows.
template <class BitmapPtr, class ValuesPtr, typename PairVector = std::vector<std::pair<uint64_t, uint8_t>>, bool ordered = true>
void call_rows(const std::vector<BitmapPtr> &columns,
               const std::vector<ValuesPtr> &column_values,
               const std::function<void(uint64_t, const PairVector &, size_t)> &callback,
               uint64_t num_rows = 0,
               bool show_progress = mtg::common::get_verbose()) {
    if (!columns.size()) {
        ProgressBar progress_bar(num_rows, "Constructed rows", std::cerr, !show_progress);
        PairVector empty;
        for (uint64_t i = 0; i < num_rows; ++i) {
            callback(i, empty, 0);
        }
        return;
    }

    if (column_values.empty())
        mtg::common::logger->warn("No column values provided, using default value of 1");

    num_rows = columns[0]->size();
    std::vector<PairVector> rows;

    ProgressBar progress_bar(num_rows, "Constructed rows", std::cerr, !show_progress);

    #pragma omp parallel for ordered num_threads(get_num_threads()) schedule(dynamic) private(rows)
    for (uint64_t i = 0; i < num_rows; i += kNumRowsInBlock) {
        uint64_t begin = i;
        uint64_t end = std::min(i + kNumRowsInBlock, num_rows);
        assert(begin <= end);

        rows.resize(end - begin);

        for (size_t j = 0; j < columns.size(); ++j) {
            if (column_values.empty()) {
                columns[j]->call_ones_in_range(begin, end, [&](uint64_t i) {
                    rows[i - begin].emplace_back(j, 1);
                });
            } else {
                using Container = std::remove_const_t<typename ValuesPtr::element_type>;
                const Container &col = *column_values[j];

                if constexpr(std::is_same_v<Container, sdsl::int_vector_buffer<>>) {
                    sdsl::int_vector_buffer<> cur_col(col.filename(), std::ios::in, col.buffersize());
                    columns[j]->call_ones_in_range(begin, end, [&](uint64_t i) {
                        rows[i - begin].emplace_back(j, cur_col[columns[j]->rank1(i) - 1]);
                    });
                } else {
                    columns[j]->call_ones_in_range(begin, end, [&](uint64_t i) {
                        rows[i - begin].emplace_back(j, col[columns[j]->rank1(i) - 1]);
                    });
                }
            }
        }

        if constexpr(ordered) {
            #pragma omp ordered
            {
                for (uint64_t j = 0; j < rows.size(); ++j) {
                    callback(j + begin, rows[j], begin / kNumRowsInBlock);
                    ++progress_bar;
                    rows[j].resize(0);
                }
            }
        } else {
            for (uint64_t j = 0; j < rows.size(); ++j) {
                callback(j + begin, rows[j], begin / kNumRowsInBlock);
                ++progress_bar;
                rows[j].resize(0);
            }
        }

    }
}

} // namespace utils

#endif // __BIT_VECTOR_UTILS_HPP__
