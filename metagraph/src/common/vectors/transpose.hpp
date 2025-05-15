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

    using Container = std::remove_const_t<typename ValuesPtr::element_type>;

    if (column_values.empty())
        mtg::common::logger->warn("No column values provided, using default value of 1");

    num_rows = columns[0]->size();

    ProgressBar progress_bar(num_rows, "Constructed rows", std::cerr, !show_progress);

    if constexpr(ordered) {
        std::vector<PairVector> rows;
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

            #pragma omp ordered
            {
                for (uint64_t j = 0; j < rows.size(); ++j) {
                    callback(j + begin, rows[j], begin / kNumRowsInBlock);
                    ++progress_bar;
                    rows[j].resize(0);
                }
            }
        }
    } else {
        size_t block_size = (num_rows + get_num_threads() - 1) / get_num_threads();
        size_t rows_per_update = 10000;
        #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
        for (uint64_t i = 0; i < num_rows; i += block_size) {
            uint64_t begin = i;
            uint64_t end = std::min(i + block_size, num_rows);
            assert(begin <= end);

            std::vector<std::shared_ptr<const Container>> cur_cols;
            cur_cols.reserve(column_values.size());
            for (const auto &col : column_values) {
                if constexpr(std::is_same_v<Container, sdsl::int_vector_buffer<>>) {
                    cur_cols.emplace_back(std::make_shared<const Container>(col->filename(),
                                                                            std::ios::in,
                                                                            col->buffersize()));
                } else {
                    cur_cols.emplace_back(std::shared_ptr<const Container>{}, col.get());
                }
            }

            std::vector<size_t> ranks(columns.size());
            std::vector<uint64_t> next_set_bit(columns.size());
            if (begin) {
                for (size_t j = 0; j < columns.size(); ++j) {
                    ranks[j] = columns[j]->rank1(begin - 1);
                }
            }

            for (size_t j = 0; j < columns.size(); ++j) {
                if (ranks[j] < columns[j]->num_set_bits()) {
                    next_set_bit[j] = columns[j]->select1(++ranks[j]);
                } else {
                    next_set_bit[j] = end;
                }
            }

            PairVector row;
            row.reserve(columns.size());
            for (size_t k = begin; k < end; ++k) {
                for (size_t j = 0; j < columns.size(); ++j) {
                    if (next_set_bit[j] == k) {
                        uint64_t val = column_values.size()
                            ? (*cur_cols[j])[ranks[j] - 1]
                            : 1;

                        row.emplace_back(j, val);
                        if (ranks[j] < columns[j]->num_set_bits()) {
                            next_set_bit[j] = columns[j]->select1(++ranks[j]);
                        } else {
                            next_set_bit[j] = end;
                        }
                    }
                }

                callback(k, row, i / block_size);
                if (k && (k % rows_per_update) == 0)
                    progress_bar += rows_per_update;

                row.clear();
            }
        }

        progress_bar += num_rows % rows_per_update;
    }
}

} // namespace utils

#endif // __BIT_VECTOR_UTILS_HPP__
