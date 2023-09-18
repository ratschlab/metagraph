#ifndef __BIT_VECTOR_UTILS_HPP__
#define __BIT_VECTOR_UTILS_HPP__

#include <vector>
#include <memory>

#include <progress_bar.hpp>

#include "bit_vector_sdsl.hpp"


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

} // namespace utils

#endif // __BIT_VECTOR_UTILS_HPP__
