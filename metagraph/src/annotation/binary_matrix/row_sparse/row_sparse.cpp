#include "row_sparse.hpp"

#include "common/logger.hpp"
#include "common/vector.hpp"
#include "graph/representation/succinct/boss.hpp"

namespace mtg {
namespace annot {
namespace binmat {

RowSparse::RowSparse(const std::function<void(const RowCallback &)> &call_rows,
                     uint64_t num_columns,
                     uint64_t num_rows,
                     uint64_t num_relations)
    : num_columns_(num_columns), num_rows_(num_columns > 0 ? num_rows : 0) {
    //TODO(ddanciu): use an int_vector_buffer instead to save memory
    // Why can't sdsl::enc_vector<> convert a 32-bit Vector??
    Vector<uint64_t> elements(num_relations);
    sdsl::bit_vector boundary(num_relations + num_rows, 0);
    uint64_t idx = 0;
    uint64_t b_idx = 0;
    call_rows([&](const auto &column_indices) {
        for (const uint64_t col : column_indices) {
            elements[idx++] = col;
            b_idx++;
        }
        boundary[b_idx++] = 1;
    });
    assert(idx == num_relations);

    boundary_ = bit_vector_small(boundary);
    elements_ = sdsl::enc_vector<>(elements);
}

bool RowSparse::get(Row row, Column column) const {
    SetBitPositions set_bits = get_row(row);
    SetBitPositions::iterator v = std::lower_bound(set_bits.begin(), set_bits.end(), column);
    return v != set_bits.end() && *v == column;
}

std::vector<BinaryMatrix::Row> RowSparse::get_column(Column column) const {
    std::vector<Row> result;
    for (Row row = 0; row < num_rows(); ++row) {
        if (get(row, column))
            result.push_back(row);
    }
    return result;
}

BinaryMatrix::SetBitPositions RowSparse::get_row(Row row) const {
    assert(boundary_[boundary_.size() - 1] == 1);
    Vector<uint64_t> result;
    uint64_t start_idx = row == 0 ? 0 : boundary_.select1(row) + 1;

    while (boundary_[start_idx] == 0) {
        result.push_back(elements_[start_idx - row]);
        start_idx++;
    }
    return result;
}

bool RowSparse::load(std::istream &f) {
    try {
        f.read(reinterpret_cast<char *>(&num_columns_), sizeof(uint64_t));
        elements_.load(f);
        boundary_.load(f);
        num_rows_ = boundary_.num_set_bits();
    } catch (...) {
        return false;
    }
    return true;
}

void RowSparse::serialize(std::ostream &f) const  {
    f.write(reinterpret_cast<const char *>(&num_columns_), sizeof(uint64_t));
    elements_.serialize(f);
    boundary_.serialize(f);
}

} // namespace binmat
} // namespace annot
} // namespace mtg
