#include "csr_matrix_flat.hpp"

#include <mutex>

#include "common/vectors/vector_algorithm.hpp"

namespace mtg::annot::matrix {

constexpr size_t buffer_size = 1024*1024;

CSRMatrixFlat::CSRMatrixFlat(const std::function<void(const std::function<void(const Vector<std::pair<uint64_t, uint64_t>>&)>&)> &row_generator,
                             size_t num_columns,
                             size_t num_rows,
                             size_t num_nonzeros,
                             uint8_t max_width) : num_columns_(num_columns) {
    sdsl::int_vector<> vals(num_nonzeros, 0, max_width);
    uint64_t r = 0;
    std::mutex mu;
    matrix_ = RowFlat<>([&](const auto &bin_row_callback) {
        row_generator([&](const auto &row) {
            std::lock_guard<std::mutex> lock(mu);
            SetBitPositions binrow;
            binrow.reserve(row.size());
            for (const auto &[j, v] : row) {
                binrow.emplace_back(j);
                vals[r++] = v;
            }

            bin_row_callback(binrow);
        });
    }, num_columns, num_rows, num_nonzeros);

    if (!max_width)
        sdsl::util::bit_compress(vals);

    values_ = std::move(vals);
}

std::vector<VectorMap<uint64_t, size_t>>
CSRMatrixFlat::get_histograms(const std::vector<size_t> &min_counts,
                              sdsl::bit_vector *unmark_discarded) const {
    bool keep_all_nonzeros = min_counts.size() == num_columns_
                                && std::all_of(min_counts.begin(), min_counts.end(), [](size_t s) { return s == 1; });
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

    for (size_t j = 0; j < num_columns_; ++j) {
        hists_maps_base[0][j][0] = n;
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
        auto &hists_map_base = hists_maps_base[thread_id];
        uint64_t flat_begin = boundaries[thread_id];
        uint64_t flat_end = boundaries[thread_id + 1];
        uint64_t begin = flat_begin / num_columns_;
        int64_t end = flat_end / num_columns_;

        auto *v_use = v_main;

        std::unique_ptr<sdsl::int_vector_buffer<>> v_copy = std::make_unique<sdsl::int_vector_buffer<>>(v_use->filename(), std::ios::in, v_use->buffersize(), 0, false, offset_);
        v_use = v_copy.get();

        uint64_t r = flat.rank1(flat_begin) - flat[flat_begin];
        auto it = v_use->begin() + r;

        uint64_t num_ones = flat_end ? flat.rank1(flat_end - 1) : 0;

        if (keep_all_nonzeros) {
            int64_t last_i = static_cast<int64_t>(begin) - 1;
            for (uint64_t r = flat_begin ? flat.rank1(flat_begin - 1) + 1 : 1; r <= num_ones; ++r) {
                uint64_t b = flat.select1(r);
                int64_t i = b / num_columns_;
                uint64_t j = b % num_columns_;
                uint64_t val = *it;
                ++it;

                if (val < HIST_CUTOFF) {
                    ++hists_map_base[j][val];
                } else {
                    ++hists_map[j][val];
                }
                ++num_nonzeros[thread_id][j];

                if (i - last_i > 1) {
                    ++last_i;
                    uint64_t num_added = i - last_i;
                    num_empty_rows.fetch_add(num_added, std::memory_order_relaxed);
                    if (unmark_discarded) {
                        for ( ; last_i < i; ++last_i) {
                            unset_bit(unmark_discarded->data(), last_i, false, std::memory_order_relaxed);
                        }
                    }
                }

                last_i = i;
            }

            if (end - last_i > 1) {
                ++last_i;
                uint64_t num_added = end - last_i;
                num_empty_rows.fetch_add(num_added, std::memory_order_relaxed);
                if (unmark_discarded) {
                    for ( ; last_i < end; ++last_i) {
                        unset_bit(unmark_discarded->data(), last_i, false, std::memory_order_relaxed);
                    }
                }
            }
        } else {
            int64_t last_i = begin;
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

            for (uint64_t r = flat_begin ? flat.rank1(flat_begin - 1) + 1 : 1; r <= num_ones; ++r) {
                uint64_t b = flat.select1(r);
                int64_t i = b / num_columns_;
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
            }
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

auto CSRMatrixFlat::get_row_values(const std::vector<BinaryMatrix::Row> &rows) const -> std::vector<RowValues> {
    std::vector<RowValues> row_values(rows.size());
    const auto &flat = matrix_.data();

    std::visit([&](const auto &v) {
        for (size_t i = 0; i < rows.size(); ++i) {
            uint64_t idx_begin = rows[i] * num_columns_;
            uint64_t idx_end = idx_begin + num_columns_;
            uint64_t r = flat.rank1(idx_begin) - flat[idx_begin];
            flat.call_ones_in_range(idx_begin, idx_end, [&](uint64_t j) {
                row_values[i].emplace_back(j % num_columns_, v[r++]);
            });
        }
    }, values_);

    return row_values;
}

bool CSRMatrixFlat::load(const std::string &fname, std::streampos offset) {
    std::ifstream fin(fname, std::ios::binary);
    fin.seekg(offset, std::ios_base::beg);
    num_columns_ = load_number(fin);
    uint64_t m_index = load_number(fin);
    uint64_t v_index = load_number(fin);
    std::ignore = m_index;
    std::ignore = v_index;

    RowFlat<> binmat;
    if (!binmat.load(fin))
        return false;

    matrix_ = std::move(binmat);
    offset_ = fin.tellg();
    values_ = sdsl::int_vector_buffer<>(fname, std::ios::in, buffer_size, 0, false, offset_);
    return true;
}

bool CSRMatrixFlat::load(std::istream &in) {
    num_columns_ = load_number(in);
    uint64_t m_index = load_number(in);
    uint64_t v_index = load_number(in);
    std::ignore = m_index;
    std::ignore = v_index;

    RowFlat<> binmat;
    if (!binmat.load(in))
        return false;

    matrix_ = std::move(binmat);
    sdsl::int_vector<> v_load;
    v_load.load(in);
    values_ = v_load;
    return true;
}

void CSRMatrixFlat::serialize(std::ostream &out) const {
    serialize_number(out, num_columns_);
    serialize_number(out, 0);
    serialize_number(out, values_.index());
    get_binary_matrix().serialize(out);
    std::visit([&](const auto &v) {
        using T = std::decay_t<decltype(v)>;
        if constexpr(std::is_same_v<T, sdsl::int_vector<>>) {
            v.serialize(out);
        }
    }, values_);
}

} // namespace mtg::annot::matrix