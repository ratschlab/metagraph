#include "rainbow.hpp"

#include <algorithm>
#include <functional>

#include <ips4o.hpp>

#include "common/serialization.hpp"
#include "common/threads/threading.hpp"
#include "common/utils/template_utils.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt.hpp"

const size_t kRowBatchSize = 1'000'000;


template <class MatrixType>
Rainbow<MatrixType>::Rainbow(MatrixType&& reduced_matrix,
                             bit_vector_rrr<>&& row_codes,
                             bit_vector_rrr<>&& row_code_delimiters,
                             uint64_t num_relations)
      : num_relations_(num_relations),
        row_codes_(std::move(row_codes)),
        row_code_delimiters_(std::move(row_code_delimiters)),
        reduced_matrix_(std::move(reduced_matrix)) {}

template <class MatrixType>
uint64_t Rainbow<MatrixType>::num_rows() const {
    return row_code_delimiters_.rank1(row_code_delimiters_.size());
}

template <class MatrixType>
uint64_t Rainbow<MatrixType>::get_code(Row row) const {
    uint64_t begin = row ? row_code_delimiters_.select1(row) + 1 : 0;
    uint64_t width = row_code_delimiters_.select1(row + 1) + 1 - begin;
    assert(width);
    assert(begin + width <= row_codes_.size());

    return row_codes_.get_int(begin, width);
}

template <class MatrixType>
bool Rainbow<MatrixType>::get(Row row, Column column) const {
    return reduced_matrix_.get(get_code(row), column);
}

template <class MatrixType>
BinaryMatrix::SetBitPositions Rainbow<MatrixType>::get_row(Row row) const {
    return reduced_matrix_.get_row(get_code(row));
}

template <class MatrixType>
std::vector<BinaryMatrix::SetBitPositions>
Rainbow<MatrixType>::get_rows(const std::vector<Row> &rows) const {
    std::vector<Row> pointers = rows;
    auto distinct_rows = get_rows(&pointers);

    std::vector<SetBitPositions> result(rows.size());
    for (size_t i = 0; i < pointers.size(); ++i) {
        result[i] = distinct_rows[pointers[i]];
    }

    return result;
}

template <class MatrixType>
std::vector<BinaryMatrix::SetBitPositions>
Rainbow<MatrixType>::get_rows(std::vector<Row> *rows, size_t num_threads) const {
    assert(rows);

    std::vector<std::pair<uint64_t, /* code */
                          uint64_t /* row */>> row_codes(rows->size());

    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < rows->size(); ++i) {
        row_codes[i] = { get_code((*rows)[i]), i };
    }

    ips4o::parallel::sort(row_codes.begin(), row_codes.end(),
                          utils::LessFirst(), num_threads);

    std::vector<Row> unique_row_codes;
    uint64_t last_code = std::numeric_limits<uint64_t>::max();

    for (const auto &[code, i] : row_codes) {
        if (code != last_code) {
            unique_row_codes.push_back(code);
            last_code = code;
        }
        (*rows)[i] = unique_row_codes.size() - 1;
    }

    std::vector<SetBitPositions> result(unique_row_codes.size());

    uint64_t batch_size = std::max((size_t)1u,
                                   std::min(kRowBatchSize,
                                            unique_row_codes.size() / num_threads));

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (uint64_t i = 0; i < unique_row_codes.size(); i += batch_size) {
        std::vector<uint64_t> indexes(unique_row_codes.begin() + i,
                                      std::min(unique_row_codes.begin() + i + batch_size,
                                               unique_row_codes.end()));
        auto rows = reduced_matrix_.get_rows(indexes);
        for (size_t j = 0; j < rows.size(); ++j) {
            result[i + j] = std::move(rows[j]);
        }
    }

    return result;
}

template <class MatrixType>
std::vector<BinaryMatrix::Row>
Rainbow<MatrixType>::get_column(Column column) const {
    sdsl::bit_vector code_column(reduced_matrix_.num_rows(), false);
    for (uint64_t r : reduced_matrix_.get_column(column)) {
        code_column[r] = true;
    }
    std::vector<Row> row_indices;
    uint64_t rows = num_rows();
    for (uint64_t i = 0; i < rows; ++i) {
        auto code = get_code(i);
        if (code_column[code])
            row_indices.emplace_back(i);
    }
    return row_indices;
}

template <class MatrixType>
bool Rainbow<MatrixType>::load(std::istream &in) {
    try {
        num_relations_ = load_number(in);
        return row_codes_.load(in)
                && row_code_delimiters_.load(in)
                && reduced_matrix_.load(in);
    } catch (...) {
        return false;
    }
}

template <class MatrixType>
void Rainbow<MatrixType>::serialize(std::ostream &out) const {
    serialize_number(out, num_relations_);
    row_codes_.serialize(out);
    row_code_delimiters_.serialize(out);
    reduced_matrix_.serialize(out);
}

template class Rainbow<BRWT>;
