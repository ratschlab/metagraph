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

auto CSRMatrixFlat::get_row_values(const std::vector<BinaryMatrix::Row> &rows) const -> std::vector<RowValues> {
    std::vector<RowValues> row_values(rows.size());
    std::visit([&](const auto &binmat) {
        using T = std::decay_t<decltype(binmat)>;
        static_assert(std::is_same_v<T, RowFlat<>>);
        const auto &flat = binmat.data();

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
    }, matrix_);
    return row_values;
}

bool CSRMatrixFlat::load(const std::string &fname, std::streampos offset) {
    std::ifstream fin(fname, std::ios::binary);
    fin.seekg(offset, std::ios_base::beg);
    num_columns_ = load_number(fin);
    uint64_t m_index = load_number(fin);
    uint64_t v_index = load_number(fin);
    std::ignore = v_index;

    switch (m_index) {
        case 0: {
            RowFlat<> binmat;
            if (!binmat.load(fin))
                return false;

            matrix_ = std::move(binmat);
        } break;
        default: {
            throw std::runtime_error("No other representations supported");
        } break;
    }

    offset_ = fin.tellg();
    values_ = sdsl::int_vector_buffer<>(fname, std::ios::in, buffer_size, 0, false, offset_);
    return true;
}

bool CSRMatrixFlat::load(std::istream &in) {
    num_columns_ = load_number(in);
    uint64_t m_index = load_number(in);
    uint64_t v_index = load_number(in);
    std::ignore = v_index;

    switch (m_index) {
        case 0: {
            RowFlat<> binmat;
            if (!binmat.load(in))
                return false;

            matrix_ = std::move(binmat);
        } break;
        default: {
            throw std::runtime_error("No other representations supported");
        } break;
    }

    sdsl::int_vector<> v_load;
    v_load.load(in);
    values_ = v_load;
    return true;
}

void CSRMatrixFlat::serialize(std::ostream &out) const {
    serialize_number(out, num_columns_);
    serialize_number(out, matrix_.index());
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