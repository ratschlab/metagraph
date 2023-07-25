#include "rainbow.hpp"

#include <algorithm>
#include <functional>

#include <ips4o.hpp>

#include "common/serialization.hpp"
#include "common/threads/threading.hpp"
#include "common/utils/template_utils.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt.hpp"


namespace mtg {
namespace annot {
namespace matrix {

template <class MatrixType>
Rainbow<MatrixType>::Rainbow(MatrixType&& reduced_matrix,
                             sdsl::bit_vector&& row_codes,
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
        row_codes_.load(in);
        return row_code_delimiters_.load(in)
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

} // namespace matrix
} // namespace annot
} // namespace mtg
