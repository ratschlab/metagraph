#include "rainbowfish.hpp"

#include <algorithm>
#include <functional>

#include <ips4o.hpp>
#include <tsl/hopscotch_map.h>

#include "common/vector.hpp"
#include "common/serialization.hpp"
#include "common/threads/threading.hpp"
#include "common/hash/hash.hpp"
#include "common/utils/template_utils.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"


namespace mtg {
namespace anno {
namespace binmat {

Rainbowfish::Rainbowfish(const std::function<void(RowCallback)> &call_rows,
                         uint64_t num_columns,
                         uint64_t buffer_size)
      : num_columns_(num_columns),
        buffer_size_(buffer_size) {
    assert(buffer_size_);

    using IndexVectorMap = tsl::hopscotch_map<SmallVector<uint32_t>,
                                              uint64_t,
                                              utils::VectorHash>;
    IndexVectorMap vector_counter;
    IndexVectorMap vector_coder;
    uint64_t total_code_length = 0;
    uint64_t num_relations_buffered = 0;

    auto flush = [&](const auto &callback) {
        if (!vector_counter.size())
            return;

        Vector<SmallVector<uint32_t>> vectors;
        vectors.reserve(vector_counter.size());

        Vector<std::pair<uint64_t, uint64_t>> index_counts;
        index_counts.reserve(vector_counter.size());

        for (auto it = vector_counter.begin(); it != vector_counter.end(); ++it) {
            index_counts.emplace_back(index_counts.size(), it->second);
            vectors.push_back(std::move(it->first));
        }
        vector_counter.clear();
        num_relations_buffered = 0;

        // sort in the order of decreasing counts
        ips4o::parallel::sort(index_counts.begin(), index_counts.end(),
            [](const auto &first, const auto &second) {
                return first.second > second.second;
            },
            get_num_threads()
        );

        SetBitPositions v;

        for (const auto &index_count : index_counts) {
            v.assign(vectors[index_count.first].begin(),
                     vectors[index_count.first].end());
            callback(v);

            assert(vector_coder.find(vectors[index_count.first]) == vector_coder.end());

            uint64_t code = vector_coder.size();

            vector_coder[std::move(vectors[index_count.first])] = code;

            total_code_length += (sdsl::bits::hi(code) + 1) * index_count.second;
        }
    };

    uint64_t rows = 0;
    SmallVector<uint32_t> row_indices_small;
    call_rows([&](const BinaryMatrix::SetBitPositions &row_indices) {
        ++rows;
        row_indices_small.assign(row_indices.begin(), row_indices.end());
        if (vector_coder.size()) {
            auto code_find = vector_coder.find(row_indices_small);
            if (code_find != vector_coder.end()) {
                total_code_length += sdsl::bits::hi(code_find.value()) + 1;
                return;
            }
        }

        auto counter_find = vector_counter.try_emplace(row_indices_small, 1);
        if (counter_find.second) {
            num_relations_buffered += row_indices_small.size();
        } else {
            counter_find.first.value()++;
        }

        if (vector_counter.size() == buffer_size_) {
            assert(num_relations_buffered <= num_columns_ * buffer_size_);
            reduced_matrix_.emplace_back(
                new ReducedMatrixType(flush,
                                      num_columns_,
                                      buffer_size_,
                                      num_relations_buffered)
            );
        }
    });

    uint64_t remaining_rows = vector_counter.size();
    assert(num_relations_buffered <= num_columns_ * remaining_rows);
    reduced_matrix_.emplace_back(
        new ReducedMatrixType(flush,
                              num_columns_,
                              remaining_rows,
                              num_relations_buffered)
    );

    sdsl::bit_vector code_bv(total_code_length);
    sdsl::bit_vector boundary_bv(total_code_length, false);

    uint64_t pos = 0;
    row_indices_small.clear();

    call_rows([&](const auto &row_indices) {
        row_indices_small.assign(row_indices.begin(), row_indices.end());

        num_relations_ += row_indices_small.size();

        uint64_t code = vector_coder[row_indices_small];
        uint64_t code_length = sdsl::bits::hi(code) + 1;

        assert(code % buffer_size_
                < reduced_matrix_[code / buffer_size_]->num_rows());

        code_bv.set_int(pos, code, code_length);
        boundary_bv[pos + code_length - 1] = true;

        pos += code_length;
    });
    vector_coder.clear();

    assert(pos == code_bv.size());

    row_code_delimiters_ = bit_vector_rrr<>(std::move(boundary_bv));
    row_codes_ = bit_vector_rrr<>(std::move(code_bv));
}

uint64_t Rainbowfish::num_rows() const {
    return row_code_delimiters_.rank1(row_code_delimiters_.size());
}

uint64_t Rainbowfish::num_distinct_rows() const {
    uint64_t num_rows = 0;

    for (const auto &block : reduced_matrix_) {
        num_rows += block->num_rows();
    }
    return num_rows;
}

uint64_t Rainbowfish::get_code(Row row) const {
    uint64_t begin = row ? row_code_delimiters_.select1(row) + 1 : 0;
    uint64_t width = row_code_delimiters_.select1(row + 1) + 1 - begin;
    assert(width);
    assert(begin + width <= row_codes_.size());

    return row_codes_.get_int(begin, width);
}

// row is in [0, num_rows), column is in [0, num_columns)
bool Rainbowfish::get(Row row, Column column) const {
    uint64_t code = get_code(row);
    return reduced_matrix_[code / buffer_size_]->get(code % buffer_size_,
                                                     column);
}

Rainbowfish::SetBitPositions Rainbowfish::get_row(Row row) const {
    uint64_t code = get_code(row);
    return reduced_matrix_[code / buffer_size_]->get_row(code % buffer_size_);
}

std::vector<Rainbowfish::SetBitPositions>
Rainbowfish::get_rows(const std::vector<Row> &rows) const {
    std::vector<std::pair<uint64_t, /* code */
                          uint64_t /* row */>> row_codes(rows.size());

    for (size_t i = 0; i < rows.size(); ++i) {
        row_codes[i] = { get_code(rows[i]), i };
    }

    std::sort(row_codes.begin(), row_codes.end(), utils::LessFirst());

    std::vector<SetBitPositions> result(rows.size());

    SetBitPositions *last_row = nullptr;
    uint64_t last_code = std::numeric_limits<uint64_t>::max();

    for (const auto &[code, row] : row_codes) {
        if (code == last_code) {
            result[row] = *last_row;
        } else {
            result[row] = reduced_matrix_[code / buffer_size_]->get_row(code % buffer_size_);
            last_row = &result[row];
            last_code = code;
        }
    }

    return result;
}

std::vector<Rainbowfish::SetBitPositions>
Rainbowfish::get_rows(std::vector<Row> *rows, size_t num_threads) const {
    assert(rows);

    std::vector<std::pair<uint64_t, /* code */
                          uint64_t /* row */>> row_codes(rows->size());

    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < rows->size(); ++i) {
        row_codes[i] = { get_code((*rows)[i]), i };
    }

    ips4o::parallel::sort(row_codes.begin(), row_codes.end(),
                          utils::LessFirst(), num_threads);

    std::vector<SetBitPositions> unique_rows;

    uint64_t last_code = std::numeric_limits<uint64_t>::max();

    // TODO: use multiple threads
    for (const auto &[code, i] : row_codes) {
        if (code != last_code) {
            unique_rows.push_back(
                reduced_matrix_[code / buffer_size_]->get_row(code % buffer_size_)
            );
            last_code = code;
        }
        (*rows)[i] = unique_rows.size() - 1;
    }

    return unique_rows;
}

std::vector<Rainbowfish::Row> Rainbowfish::get_column(Column column) const {
    sdsl::bit_vector distinct_row_indices(num_distinct_rows(), false);
    uint64_t offset = 0;
    for (const auto &mat : reduced_matrix_) {
        for (const auto &a : mat->get_column(column)) {
            distinct_row_indices[offset + a] = true;
        }
        offset += mat->num_rows();
    }
    std::vector<Row> row_indices;
    uint64_t rows = num_rows();
    for (uint64_t i = 0; i < rows; ++i) {
        auto code = get_code(i);
        if (distinct_row_indices[code])
            row_indices.emplace_back(i);
    }
    return row_indices;
}

bool Rainbowfish::load(std::istream &in) {
    try {
        num_columns_ = load_number(in);
        num_relations_ = load_number(in);
        buffer_size_ = load_number(in);
        if (!row_codes_.load(in) || !row_code_delimiters_.load(in))
            return false;

        // TODO: handle diff types
        reduced_matrix_.clear();
        uint64_t num_blocks = load_number(in);
        reduced_matrix_.reserve(num_blocks);
        while (num_blocks--) {
            reduced_matrix_.emplace_back(new ReducedMatrixType());
            reduced_matrix_.back()->load(in);
        }
    } catch (...) {
        return false;
    }

    return true;
}

void Rainbowfish::serialize(std::ostream &out) const {
    serialize_number(out, num_columns_);
    serialize_number(out, num_relations_);
    serialize_number(out, buffer_size_);
    row_codes_.serialize(out);
    row_code_delimiters_.serialize(out);

    // TODO: handle diff types
    serialize_number(out, reduced_matrix_.size());
    for (auto &matrix : reduced_matrix_) {
        matrix->serialize(out);
    }
}

// number of ones in the matrix
uint64_t Rainbowfish::num_relations() const {
    return num_relations_;
}

} // namespace binmat
} // namespace anno
} // namespace mtg
