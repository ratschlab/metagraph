#ifndef __CSR_MATRIX_FLAT_HPP__
#define __CSR_MATRIX_FLAT_HPP__

#include <variant>
#include <omp.h>

#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <tsl/hopscotch_map.h>
#include <progress_bar.hpp>

#include "annotation/int_matrix/base/int_matrix.hpp"
#include "annotation/binary_matrix/row_flat/flat_matrix.hpp"
#include "annotation/binary_matrix/row_disk/row_disk.hpp"
#include "common/serialization.hpp"


namespace mtg::annot::matrix {

class CSRMatrixFlat : public RowMajor, public IntMatrix {
  public:
    CSRMatrixFlat() {}

    CSRMatrixFlat(const std::function<void(const std::function<void(const Vector<std::pair<uint64_t, uint64_t>>&)>&)> &row_generator,
                  size_t num_columns,
                  size_t num_rows,
                  size_t num_nonzeros,
                  uint8_t max_width = 0);

    template <class BinMat, class Values>
    CSRMatrixFlat(BinMat&& binmat, Values&& values)
        : matrix_(std::move(binmat)), values_(std::move(values)) {}

    std::vector<RowValues> get_row_values(const std::vector<BinaryMatrix::Row> &rows) const;

    void call_row_values(const std::function<void(uint64_t, RowValues&&, size_t)> &callback,
                         bool ordered = true) const override {
        constexpr size_t kNumRowsInBlock = 50'000;
        std::visit([&](const auto &binmat) {
            using T = std::decay_t<decltype(binmat)>;
            static_assert(std::is_same_v<T, RowFlat<>>);
            const auto &flat = binmat.data();

            std::visit([&](const auto &v) {
                using S = std::decay_t<decltype(v)>;
                uint64_t n = num_rows();

                if (ordered) {
                    ProgressBar progress_bar(n, "Constructed rows", std::cerr, !common::get_verbose());
                    uint64_t block_size = kNumRowsInBlock - (kNumRowsInBlock % num_columns_);
                    std::vector<RowValues> rows;
                    #pragma omp parallel for ordered num_threads(get_num_threads()) schedule(dynamic) private(rows)
                    for (uint64_t begin = 0; begin < n; begin += block_size) {
                        size_t thread_id = begin / block_size;
                        std::unique_ptr<sdsl::int_vector_buffer<>> v_copy;
                        if constexpr(std::is_same_v<S, sdsl::int_vector_buffer<>>) {
                            v_copy = std::make_unique<sdsl::int_vector_buffer<>>(v.filename(), std::ios::in, v.buffersize(), 0, false, offset_);
                        }

                        uint64_t end = std::min(begin + block_size, n);
                        assert(begin <= end);
                        rows.resize(end - begin);

                        uint64_t flat_begin = begin * num_columns_;
                        uint64_t flat_end = end * num_columns_;
                        uint64_t r = flat.rank1(flat_begin) - flat[flat_begin];

                        auto set_bit_callback = [&](uint64_t b) __attribute__((always_inline)) {
                            uint64_t i = b / num_columns_;
                            uint64_t j = b % num_columns_;
                            uint64_t val = 0;

                            if constexpr(std::is_same_v<S, sdsl::int_vector<>>) {
                                val = v[r++];
                            }

                            if constexpr(std::is_same_v<S, sdsl::int_vector_buffer<>>) {
                                val = (*v_copy)[r++];
                            }

                            rows[i - begin].emplace_back(j, val);
                        };

                        flat.call_ones_in_range(flat_begin, flat_end, set_bit_callback);

                        #pragma omp ordered
                        {
                            uint64_t i = begin;
                            progress_bar += rows.size();
                            for (auto &row : rows) {
                                callback(i++, std::move(row), thread_id);
                                row = RowValues();
                            }
                        }
                    }
                } else {
                    ProgressBar progress_bar(n, "Constructed rows", std::cerr, !common::get_verbose());
                    uint64_t block_size = n / std::max(1u, get_num_threads() - 1);
                    block_size -= block_size % 64;
                    #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
                    for (uint64_t begin = 0; begin < n; begin += block_size) {
                        size_t thread_id = begin / block_size;
                        RowValues row;

                        const S *v_use = &v;

                        std::unique_ptr<sdsl::int_vector_buffer<>> v_copy;
                        if constexpr(std::is_same_v<S, sdsl::int_vector_buffer<>>) {
                            v_copy = std::make_unique<sdsl::int_vector_buffer<>>(v.filename(), std::ios::in, v.buffersize(), 0, false, offset_);
                            v_use = v_copy.get();
                        }

                        uint64_t end = std::min(begin + block_size, n);
                        assert(begin <= end);
                        uint64_t flat_begin = begin * num_columns_;
                        uint64_t flat_end = end * num_columns_;
                        uint64_t r = flat.rank1(flat_begin) - flat[flat_begin];
                        auto it = const_cast<S*>(v_use)->begin() + r;

                        uint64_t last_i = begin;
                        auto set_bit_callback = [&](uint64_t b) __attribute__((always_inline)) {
                            uint64_t i = b / num_columns_;
                            if (i > last_i) {
                                callback(last_i, std::move(row), thread_id);

                                ++last_i;
                                row = RowValues();

                                for ( ; last_i < i; ++last_i) {
                                    callback(last_i, RowValues(), thread_id);
                                }
                            }

                            uint64_t j = b % num_columns_;
                            uint64_t val = *it;
                            ++it;

                            row.emplace_back(j, val);
                        };

                        flat.call_ones_in_range(flat_begin, flat_end, set_bit_callback);

                        callback(last_i, std::move(row), thread_id);
                        ++last_i;

                        for ( ; last_i < end; ++last_i) {
                            callback(last_i, RowValues(), thread_id);
                        }

                        progress_bar += end - begin;
                    }
                }
            }, values_);
        }, matrix_);
    }

    std::vector<VectorMap<uint64_t, size_t>> get_histograms(const std::vector<size_t> &min_counts = {},
                                                            sdsl::bit_vector *unmark_discarded = nullptr) const override {
        constexpr uint8_t HIST_CUTOFF = 64;
        // bool parallel = get_num_threads() > 1;
        uint64_t n = num_rows();
        uint64_t block_size = n / std::max(1u, get_num_threads() - 1);
        block_size -= block_size % 64;
        uint64_t num_tasks = (n + block_size - 1) / block_size;
        std::vector<std::vector<std::vector<size_t>>> hists_maps_base;
        std::vector<std::vector<tsl::hopscotch_map<uint64_t, size_t>>> hists_maps;
        for (size_t j = 0; j < num_tasks; ++j) {
            uint64_t cur_block_size = block_size;
            if ((j + 1) * block_size > n)
                cur_block_size -= (j + 1) * block_size - n;

            hists_maps.emplace_back(get_binary_matrix().num_columns());
            hists_maps_base.emplace_back(get_binary_matrix().num_columns());
            for (size_t i = 0; i < num_columns_; ++i) {
                hists_maps_base.back()[i].resize(HIST_CUTOFF);
                hists_maps_base.back()[i][0] = cur_block_size;
            }
        }

        std::atomic<uint64_t> num_empty_rows{0};

        std::atomic_thread_fence(std::memory_order_release);

        auto row_callback = [&](uint64_t row_i, auto&& row, size_t bucket_id) __attribute__((always_inline)) {
            bool keep = min_counts.empty();

            if (min_counts.size()) {
                for (const auto &[j, c] : row) {
                    if (c >= min_counts[j]) {
                        keep = true;
                        break;
                    }
                }
            }

            if (!keep) {
                // TODO: are these bottlenecks?
                if (unmark_discarded)
                    unset_bit(unmark_discarded->data(), row_i, false);
                    // unset_bit(unmark_discarded->data(), row_i, parallel, std::memory_order_relaxed);

                num_empty_rows.fetch_add(1, std::memory_order_relaxed);
            } else {
                auto &hists_map = hists_maps[bucket_id];
                auto &hists_map_base = hists_maps_base[bucket_id];
                for (const auto &[j, c] : row) {
                    if (c < HIST_CUTOFF) {
                        ++hists_map_base[j][c];
                    } else {
                        ++hists_map[j][c];
                    }
                    --hists_map_base[j][0];
                }
            }
        };

        call_row_values(row_callback, false);

        std::atomic_thread_fence(std::memory_order_acquire);

        std::vector<VectorMap<uint64_t, size_t>> hist_maps(num_columns_);
        common::logger->trace("Merging histograms");
        for (size_t j = 0; j < num_columns_; ++j) {
            for (const auto &[k, c] : hists_maps[0][j]) {
                hist_maps[j][k] += c;
            }

            for (size_t i = 1; i < hists_maps_base.size(); ++i) {
                for (size_t k = 0; k < HIST_CUTOFF; ++k) {
                    hists_maps_base[0][j][k] += hists_maps_base[i][j][k];
                }

                for (const auto &[k, c] : hists_maps[i][j]) {
                    hist_maps[j][k] += c;
                }
            }
        }

        for (size_t j = 0; j < num_columns_; ++j) {
            hists_maps_base[0][j][0] -= num_empty_rows;
            for (size_t k = 0; k < HIST_CUTOFF; ++k) {
                if (hists_maps_base[0][j][k] > 0)
                    hist_maps[j][k] += hists_maps_base[0][j][k];
            }
        }

        return hist_maps;
    }

    uint64_t num_columns() const { return num_columns_; }
    uint64_t num_rows() const { return get_row_major_binary_matrix().num_rows(); }
    uint64_t num_relations() const { return get_row_major_binary_matrix().num_relations(); }
    SetBitPositions get_row(Row row) const { return get_row_major_binary_matrix().get_row(row); }
    std::vector<Row> get_column(Column column) const { return get_row_major_binary_matrix().get_column(column); }

    const BinaryMatrix& get_binary_matrix() const { return get_row_major_binary_matrix(); }

    bool load(std::istream &in);
    bool load(const std::string &fname, std::streampos offset = 0);
    void serialize(std::ostream &out) const;

  private:
    size_t num_columns_ = 0;

    // using BinMatStorage = std::variant<RowFlat<>, RowDisk>;
    using BinMatStorage = std::variant<RowFlat<>>;
    using ValueStorage = std::variant<sdsl::int_vector<>, sdsl::int_vector_buffer<>>;

    BinMatStorage matrix_;
    ValueStorage values_;

    std::streampos offset_;

    const RowMajor& get_row_major_binary_matrix() const {
        return std::visit([&](const auto &m) -> const RowMajor& { return m; }, matrix_);
    }
};

} // namespace mtg::annot::matrix

#endif // __CSR_MATRIX_FLAT_HPP__