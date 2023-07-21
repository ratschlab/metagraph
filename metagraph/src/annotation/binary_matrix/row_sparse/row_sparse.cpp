#include "row_sparse.hpp"

#include "common/logger.hpp"
#include "common/utils/file_utils.hpp"
#include "graph/representation/succinct/boss.hpp"

namespace mtg {
namespace annot {
namespace matrix {

RowSparse::RowSparse(const std::function<void(const RowCallback &)> &call_rows,
                     uint64_t num_columns,
                     uint64_t num_rows,
                     uint64_t num_relations)
      : num_columns_(num_columns),
        num_rows_(num_columns > 0 ? num_rows : 0) {
    uint8_t col_index_width = num_columns ? sdsl::bits::hi(num_columns - 1) + 1 : 1;
    const auto tmp_dir = utils::create_temp_dir(utils::get_swap_path(), "set_bits");
    sdsl::int_vector_buffer<> set_bits(tmp_dir/"vector", std::ios::out,
                                       1024 * 1024, col_index_width);
    sdsl::bit_vector boundary(num_relations + num_rows, 0);
    uint64_t b_idx = 0;
    call_rows([&](const auto &column_indices) {
        assert(std::is_sorted(column_indices.begin(), column_indices.end()));
        uint64_t last = 0;
        for (uint64_t col : column_indices) {
            set_bits.push_back(col - last);
            last = col;
            b_idx++;
        }
        boundary[b_idx++] = 1;
    });
    assert(set_bits.size() == num_relations);

    boundary_ = bit_vector_small(boundary);
    set_bits_ = sdsl::vlc_vector<>(set_bits);
    utils::remove_temp_dir(tmp_dir);
}

std::vector<BinaryMatrix::Row> RowSparse::get_column(Column column) const {
    std::vector<Row> result;
    for (Row row = 0; row < num_rows(); ++row) {
        SetBitPositions set_bits = get_row(row);
        if (std::binary_search(set_bits.begin(), set_bits.end(), column))
            result.push_back(row);
    }
    return result;
}

BinaryMatrix::SetBitPositions RowSparse::get_row(Row row) const {
    assert(boundary_[boundary_.size() - 1] == 1);
    uint64_t start_idx = row == 0 ? 0 : boundary_.select1(row) + 1;
    uint64_t end_idx = boundary_.next1(start_idx);
    SetBitPositions result;
    result.reserve(end_idx - start_idx);
    // In each row, the first value in `set_bits_` stores the first set bit,
    // and all next ones store deltas pos[i] - pos[i-1].
    uint64_t last = 0;
    for (uint64_t i = start_idx; i != end_idx; ++i) {
        result.push_back(last += set_bits_[i - row]);
    }
    return result;
}

bool RowSparse::load(std::istream &f) {
    try {
        f.read(reinterpret_cast<char *>(&num_columns_), sizeof(uint64_t));
        set_bits_.load(f);
        boundary_.load(f);
        num_rows_ = boundary_.num_set_bits();
    } catch (...) {
        return false;
    }
    return true;
}

void RowSparse::serialize(std::ostream &f) const  {
    f.write(reinterpret_cast<const char *>(&num_columns_), sizeof(uint64_t));
    set_bits_.serialize(f);
    boundary_.serialize(f);
}

} // namespace matrix
} // namespace annot
} // namespace mtg
