#include "binary_matrix.hpp"

#include <ips4o.hpp>
#include <tsl/hopscotch_map.h>

#include "common/vectors/bitmap.hpp"
#include "common/serialization.hpp"
#include "common/utils/template_utils.hpp"
#include "common/vector_map.hpp"


namespace mtg {
namespace annot {
namespace binmat {

std::vector<BinaryMatrix::SetBitPositions>
BinaryMatrix::get_rows(const std::vector<Row> &row_ids) const {
    std::vector<SetBitPositions> rows(row_ids.size());

    auto slice = slice_rows(row_ids);

    assert(slice.size() >= row_ids.size());

    auto row_begin = slice.begin();

    for (size_t i = 0; i < rows.size(); ++i) {
        // every row in `slice` ends with `-1`
        auto row_end = std::find(row_begin, slice.end(),
                                 std::numeric_limits<Column>::max());
        rows[i].assign(row_begin, row_end);
        row_begin = row_end + 1;
    }

    return rows;
}

void BinaryMatrix::call_columns(const std::vector<Column> &column_ids,
                                const std::function<void(size_t, const bitmap&)> &callback,
                                size_t num_threads) const {
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (size_t k = 0; k < column_ids.size(); ++k) {
        callback(k, bitmap_generator(get_column(column_ids[k]), num_rows()));
    }
}

std::vector<std::pair<BinaryMatrix::Column, size_t /* count */>>
BinaryMatrix::sum_rows(const std::vector<std::pair<Row, size_t>> &index_counts,
                       size_t min_count) const {
    min_count = std::max(min_count, size_t(1));

    size_t total_sum_count = 0;
    for (auto [i, count] : index_counts) {
        total_sum_count += count;
    }

    if (total_sum_count < min_count)
        return {};

    // TODO: call slice_rows

    VectorMap<Column, size_t> col_counts;

    for (auto [i, count] : index_counts) {
        for (size_t j : get_row(i)) {
            col_counts[j] += count;
        }
    }

    auto &result = const_cast<std::vector<std::pair<Column, size_t>>&>(col_counts.values_container());

    result.erase(std::remove_if(result.begin(), result.end(),
                                [&](const auto &p) { return p.second < min_count; }),
                 result.end());

    return std::move(result);
}


std::vector<RainbowMatrix::SetBitPositions>
RainbowMatrix::get_rows(std::vector<Row> *rows, size_t num_threads) const {
    assert(rows);

    std::vector<std::pair<uint64_t, /* code */
                          uint64_t /* row */>> row_codes(rows->size());

    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < rows->size(); ++i) {
        row_codes[i] = { get_code((*rows)[i]), i };
    }

    ips4o::parallel::sort(row_codes.begin(), row_codes.end(),
                          utils::LessFirst(), num_threads);

    std::vector<uint64_t> codes;
    codes.reserve(row_codes.size());
    uint64_t last_code = std::numeric_limits<uint64_t>::max();
    for (const auto &[code, i] : row_codes) {
        if (code != last_code) {
            codes.push_back(code);
            last_code = code;
        }
        (*rows)[i] = codes.size() - 1;
    }
    row_codes = {};

    std::vector<SetBitPositions> unique_rows(codes.size());

    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < codes.size(); ++i) {
        unique_rows[i] = code_to_row(codes[i]);
    }

    return unique_rows;
}

// TODO: improve
RainbowMatrix::SetBitPositions
RainbowMatrix::slice_rows(const std::vector<Row> &rows) const {
    SetBitPositions slice;
    slice.reserve(rows.size() * 2);

    for (const auto &row : get_rows(rows)) {
        for (uint64_t j : row) {
            slice.push_back(j);
        }
        slice.push_back(std::numeric_limits<Column>::max());
    }

    return slice;
}

std::vector<RainbowMatrix::SetBitPositions>
RainbowMatrix::get_rows(const std::vector<Row> &rows) const {
    std::vector<Row> pointers = rows;
    auto distinct_rows = get_rows(&pointers);

    std::vector<SetBitPositions> result(rows.size());
    for (size_t i = 0; i < pointers.size(); ++i) {
        result[i] = distinct_rows[pointers[i]];
    }

    return result;
}

// TODO: merge with BinaryMatrix::sum_rows
std::vector<std::pair<RainbowMatrix::Column, size_t /* count */>>
RainbowMatrix::sum_rows(const std::vector<std::pair<Row, size_t>> &index_counts,
                        size_t min_count) const {
    min_count = std::max(min_count, size_t(1));

    size_t total_sum_count = 0;
    for (const auto &pair : index_counts) {
        total_sum_count += pair.second;
    }

    if (total_sum_count < min_count)
        return {};

    tsl::hopscotch_map<size_t, size_t> code_count;
    code_count.reserve(index_counts.size());
    for (auto [i, count] : index_counts) {
        code_count[get_code(i)] += count;
    }

    VectorMap<Column, size_t> col_counts;

    for (auto [c, count] : code_count) {
        for (size_t j : code_to_row(c)) {
            col_counts[j] += count;
        }
    }

    auto &result = const_cast<std::vector<std::pair<Column, size_t>>&>(col_counts.values_container());

    result.erase(std::remove_if(result.begin(), result.end(),
                                [&](const auto &p) { return p.second < min_count; }),
                 result.end());

    return std::move(result);
}


template <typename RowType>
StreamRows<RowType>::StreamRows(const std::string &filename, size_t offset) {
    std::ifstream instream(filename, std::ios::binary);

    if (!instream.good() || !instream.seekg(offset).good())
        throw std::ifstream::failure("Cannot read rows from file " + filename);

    (void)load_number(instream);
    (void)load_number(instream);

    inbuf_ = sdsl::int_vector_buffer<>(filename,
                                       std::ios::in | std::ios::binary,
                                       1024 * 1024,
                                       0,
                                       false,
                                       instream.tellg());
}

template <typename RowType>
RowType* StreamRows<RowType>::next_row() {
    row_.clear();

    while (i_ < inbuf_.size()) {
        auto value = inbuf_[i_++];
        if (value - 1 > std::numeric_limits<typename RowType::value_type>::max())
            throw std::ifstream::failure("Integer overflow: trying to read too"
                                         " large column index: " + std::to_string(value - 1));
        if (value) {
            row_.push_back(value - 1);
        } else {
            return &row_;
        }
    }
    return nullptr;
}

template class StreamRows<BinaryMatrix::SetBitPositions>;


void append_row_major(const std::string &filename,
                      const std::function<void(BinaryMatrix::RowCallback)> &call_rows,
                      uint64_t num_cols) {
    std::ofstream outstream(filename, std::ios::binary | std::ios::app);

    uint64_t num_rows = 0;

    // write dummy num_rows value to fill in later
    const uint64_t header_offs = outstream.tellp();
    serialize_number(outstream, 0);
    serialize_number(outstream, num_cols);
    const uint64_t iv_offs = outstream.tellp();
    outstream.close();

    {
        auto outbuf = sdsl::int_vector_buffer<>(filename,
                                                std::ios::out | std::ios::binary,
                                                1024 * 1024,
                                                sdsl::bits::hi(num_cols) + 1,
                                                false,
                                                iv_offs);

        call_rows([&](const auto &row) {
            for (auto val : row) {
                outbuf.push_back(val + 1);
            }
            outbuf.push_back(0);
            num_rows++;
        });

        outbuf.close();
    }

    outstream.open(filename, std::ios::in | std::ios::out | std::ios::binary);
    outstream.seekp(header_offs);
    serialize_number(outstream, num_rows);
}

} // namespace binmat
} // namespace annot
} // namespace mtg
