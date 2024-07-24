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

    inline void call_row_values(const std::function<void(uint64_t, RowValues&&, size_t)> &callback,
                                bool ordered = true) const override final {
        constexpr size_t kNumRowsInBlock = 50'000;
        const auto &flat = matrix_.data();

        sdsl::int_vector_buffer<> *v_main;
        std::visit([&](const auto &v) {
            using S = std::decay_t<decltype(v)>;
            if constexpr(std::is_same_v<S, sdsl::int_vector_buffer<>>) {
                v_main = const_cast<S*>(&v);
            } else {
                throw std::runtime_error("Failed");
            }
        }, values_);

        uint64_t n = num_rows();

        if (ordered) {
            ProgressBar progress_bar(n, "Constructed rows", std::cerr, !common::get_verbose());
            uint64_t block_size = kNumRowsInBlock - (kNumRowsInBlock % num_columns_);
            std::vector<RowValues> rows;
            #pragma omp parallel for ordered num_threads(get_num_threads()) schedule(dynamic) private(rows)
            for (uint64_t begin = 0; begin < n; begin += block_size) {
                size_t thread_id = begin / block_size;
                std::unique_ptr<sdsl::int_vector_buffer<>> v_copy = std::make_unique<sdsl::int_vector_buffer<>>(v_main->filename(), std::ios::in, v_main->buffersize(), 0, false, offset_);

                uint64_t end = std::min(begin + block_size, n);
                assert(begin <= end);
                rows.resize(end - begin);

                uint64_t flat_begin = begin * num_columns_;
                uint64_t flat_end = end * num_columns_;
                uint64_t r = flat.rank1(flat_begin) - flat[flat_begin];
                auto it = v_copy->begin() + r;

                flat.call_ones_in_range(flat_begin, flat_end, [&](uint64_t b) __attribute__((always_inline)) {
                    uint64_t i = b / num_columns_;
                    uint64_t j = b % num_columns_;
                    uint64_t val = *it;
                    ++it;
                    rows[i - begin].emplace_back(j, val);
                });

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
            uint64_t num_set_bits = flat.num_set_bits();
            uint64_t block_size = num_set_bits / std::max(1u, get_num_threads() - 1);
            size_t row_block_size = num_columns_ * 64;
            std::vector<uint64_t> boundaries;
            boundaries.emplace_back(0);
            for (size_t r = 1; r < num_set_bits; r += block_size) {
                uint64_t b = flat.select1(r);
                uint64_t b_round_ncols_64 = (b + row_block_size - 1) / row_block_size * row_block_size;
                if (b_round_ncols_64 > boundaries.back()) {
                    boundaries.emplace_back(b_round_ncols_64);
                    if (boundaries.back() >= flat.size()) {
                        boundaries.back() = flat.size();
                        break;
                    }
                }
            }

            if (boundaries.back() != flat.size())
                boundaries.emplace_back(flat.size());

            common::logger->trace("Streaming rows");
            ProgressBar progress_bar(n, "Constructed rows", std::cerr, !common::get_verbose());
            #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
            for (size_t thread_id = 0; thread_id < boundaries.size() - 1; ++thread_id) {
                uint64_t flat_begin = boundaries[thread_id];
                uint64_t flat_end = boundaries[thread_id + 1];
                uint64_t begin = flat_begin / num_columns_;
                uint64_t end = flat_end / num_columns_;

                RowValues row;

                auto *v_use = v_main;

                std::unique_ptr<sdsl::int_vector_buffer<>> v_copy = std::make_unique<sdsl::int_vector_buffer<>>(v_use->filename(), std::ios::in, v_use->buffersize(), 0, false, offset_);
                v_use = v_copy.get();

                uint64_t r = flat.rank1(flat_begin) - flat[flat_begin];
                auto it = v_use->begin() + r;

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
    }

    inline std::vector<VectorMap<uint64_t, size_t>>
    get_histograms(const std::vector<size_t> &min_counts = {},
                   sdsl::bit_vector *unmark_discarded = nullptr) const override final {
        constexpr uint8_t HIST_CUTOFF = 64;
        uint64_t n = num_rows();
        size_t num_tasks = get_num_threads() + 1;
        std::vector<std::vector<size_t>> num_nonzeros;
        std::vector<std::vector<std::vector<size_t>>> hists_maps_base;
        for (size_t j = 0; j < num_tasks; ++j) {
            hists_maps_base.emplace_back(get_binary_matrix().num_columns());
            num_nonzeros.emplace_back(get_binary_matrix().num_columns());
            for (size_t i = 0; i < num_columns_; ++i) {
                hists_maps_base.back()[i].resize(HIST_CUTOFF);
            }
        }

        for (size_t i = 0; i < num_columns_; ++i) {
            hists_maps_base[0][i][0] = n;
        }

        std::atomic<uint64_t> num_empty_rows{0};

        std::atomic_thread_fence(std::memory_order_release);

        const auto &flat = matrix_.data();

        sdsl::int_vector_buffer<> *v_main;
        std::visit([&](const auto &v) {
            using S = std::decay_t<decltype(v)>;
            if constexpr(std::is_same_v<S, sdsl::int_vector_buffer<>>) {
                v_main = const_cast<S*>(&v);
            } else {
                throw std::runtime_error("Failed");
            }
        }, values_);

        uint64_t num_set_bits = flat.num_set_bits();
        uint64_t block_size = num_set_bits / std::max(1u, get_num_threads() - 1);
        size_t row_block_size = num_columns_ * 64;
        std::vector<uint64_t> boundaries;
        boundaries.emplace_back(0);
        for (size_t r = 1; r < num_set_bits; r += block_size) {
            uint64_t b = flat.select1(r);
            uint64_t b_round_ncols_64 = (b + row_block_size - 1) / row_block_size * row_block_size;
            if (b_round_ncols_64 > boundaries.back()) {
                boundaries.emplace_back(b_round_ncols_64);
                if (boundaries.back() >= flat.size()) {
                    boundaries.back() = flat.size();
                    break;
                }
            }
        }

        if (boundaries.back() != flat.size())
            boundaries.emplace_back(flat.size());

        std::vector<VectorMap<uint64_t, size_t>> hist_maps(num_columns_);

        common::logger->trace("Streaming rows");
        ProgressBar progress_bar(n, "Constructed rows", std::cerr, !common::get_verbose());
        #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
        for (size_t thread_id = 0; thread_id < boundaries.size() - 1; ++thread_id) {
            std::vector<tsl::hopscotch_map<uint64_t, size_t>> hists_map(get_binary_matrix().num_columns());
            std::vector<std::pair<uint64_t, size_t>> row;
            auto row_callback = [&](uint64_t row_i) __attribute__((always_inline)) {
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
                    if (unmark_discarded)
                        unset_bit(unmark_discarded->data(), row_i, false, std::memory_order_relaxed);

                    num_empty_rows.fetch_add(1, std::memory_order_relaxed);
                } else {
                    auto &hists_map_base = hists_maps_base[thread_id];
                    for (const auto &[j, c] : row) {
                        if (c < HIST_CUTOFF) {
                            ++hists_map_base[j][c];
                        } else {
                            ++hists_map[j][c];
                        }
                        ++num_nonzeros[thread_id][j];
                    }
                }
            };
            uint64_t flat_begin = boundaries[thread_id];
            uint64_t flat_end = boundaries[thread_id + 1];
            uint64_t begin = flat_begin / num_columns_;
            uint64_t end = flat_end / num_columns_;

            auto *v_use = v_main;

            std::unique_ptr<sdsl::int_vector_buffer<>> v_copy = std::make_unique<sdsl::int_vector_buffer<>>(v_use->filename(), std::ios::in, v_use->buffersize(), 0, false, offset_);
            v_use = v_copy.get();

            uint64_t r = flat.rank1(flat_begin) - flat[flat_begin];
            auto it = v_use->begin() + r;

            uint64_t last_i = begin;
            uint64_t num_ones = flat_end ? flat.rank1(flat_end - 1) : 0;
            for (uint64_t r = flat_begin ? flat.rank1(flat_begin - 1) + 1 : 1; r <= num_ones; ++r) {
                uint64_t b = flat.select1(r);
                uint64_t i = b / num_columns_;
                if (i > last_i) {
                    row_callback(last_i);

                    ++last_i;
                    row.clear();

                    num_empty_rows.fetch_add(i - last_i, std::memory_order_relaxed);
                    if (unmark_discarded) {
                        for ( ; last_i < i; ++last_i) {
                            unset_bit(unmark_discarded->data(), last_i, false, std::memory_order_relaxed);
                        }
                    } else {
                        last_i = i;
                    }
                }

                uint64_t j = b % num_columns_;
                uint64_t val = *it;
                ++it;

                row.emplace_back(j, val);
            }

            row_callback(last_i);
            ++last_i;

            num_empty_rows.fetch_add(end - last_i, std::memory_order_relaxed);
            if (unmark_discarded) {
                for ( ; last_i < end; ++last_i) {
                    unset_bit(unmark_discarded->data(), last_i, false, std::memory_order_relaxed);
                }
            } else {
                last_i = end;
            }

            progress_bar += end - begin;

            #pragma omp critical
            {
                for (size_t j = 0; j < num_columns_; ++j) {
                    for (const auto &[k, c] : hists_map[j]) {
                        hist_maps[j][k] += c;
                    }
                }
            }
        }

        std::atomic_thread_fence(std::memory_order_acquire);


        common::logger->trace("Merging histograms");
        for (size_t j = 0; j < num_columns_; ++j) {
            for (size_t i = 1; i < hists_maps_base.size(); ++i) {
                for (size_t k = 0; k < HIST_CUTOFF; ++k) {
                    hists_maps_base[0][j][k] += hists_maps_base[i][j][k];
                }
            }
        }

        for (size_t j = 0; j < num_columns_; ++j) {
            hists_maps_base[0][j][0] -= num_empty_rows;
            for (size_t i = 0; i < num_nonzeros.size(); ++i) {
                hists_maps_base[0][j][0] -= num_nonzeros[i][j];
            }
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

    using BitMatStorage = RowFlat<>;
    using ValueStorage = std::variant<sdsl::int_vector<>, sdsl::int_vector_buffer<>>;

    BitMatStorage matrix_;
    ValueStorage values_;

    std::streampos offset_;

    const RowMajor& get_row_major_binary_matrix() const {
        return matrix_;
    }
};

} // namespace mtg::annot::matrix

#endif // __CSR_MATRIX_FLAT_HPP__