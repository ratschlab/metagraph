#include "bin_rel_wt_sdsl.hpp"

#include <cassert>
#include <iterator>
#include <cmath>

#include "common/serialization.hpp"
#include "common/vectors/vector_algorithm.hpp"


BinRelWT_sdsl
::BinRelWT_sdsl(const std::function<void(const RowCallback &)> &generate_rows,
                uint64_t num_relations,
                uint64_t num_columns)
      : num_columns_(num_columns) {
    sdsl::int_vector<> flat(num_relations, 0, sdsl::bits::hi(num_columns) + 1);

    // delimiters_ includes a 0 for each char in base
    // string and a 1 for each row (objects).
    // delimiter_ starts and ends with 1.
    std::vector<bool> delimiters_vec = { 1, };

    uint64_t index = 0;

    generate_rows([&](const SetBitPositions &row_set_bits) {
        for (const auto &col_index : row_set_bits) {
            assert(col_index < num_columns);
            assert(index < flat.size());
            flat[index++] = col_index;
            delimiters_vec.push_back(0);
        }
        delimiters_vec.push_back(1);
    });

    assert(index == num_relations);

    delimiters_ = bit_vector_rrr<>(to_sdsl(std::move(delimiters_vec)));

    decltype(wt_)(std::move(flat)).swap(wt_);
}

uint64_t BinRelWT_sdsl::num_columns() const {
    return num_columns_;
}

uint64_t BinRelWT_sdsl::num_rows() const {
    assert(delimiters_.num_set_bits());
    return delimiters_.num_set_bits() - 1;
}

BinRelWT_sdsl::SetBitPositions BinRelWT_sdsl::get_row(Row row) const {
    // delimiters_ starts with a 1 indicating the
    // beginning of first row and ends with a 1,
    // indicating the end of last line.
    // Note that wt_ stores a character for each 0 in delimiters_
    // and no character for each 1 in delimiters_.

    // Not counting i first 1s indicating the boundries.
    uint64_t first_string_index = delimiters_.select1(row + 1) - row;

    // Not counting i+1 first 1s indicating the boundries.
    uint64_t last_string_index = delimiters_.select1(row + 2) - (row + 1);

    // Get label indices from the base string stored in wt_.
    typedef sdsl::int_vector<>::size_type size_type;
    typedef typename decltype(wt_)::value_type value_type;

    size_type num_row_set_bits;

    std::vector<size_type> rank_c_i(last_string_index - first_string_index);
    std::vector<size_type> rank_c_j(last_string_index - first_string_index);

    std::vector<value_type> label_indices(last_string_index - first_string_index);

    wt_.interval_symbols(first_string_index,
                         last_string_index,
                         num_row_set_bits,
                         label_indices,
                         rank_c_i,
                         rank_c_j);

    return SetBitPositions(label_indices.begin(), label_indices.end());
}

std::vector<BinRelWT_sdsl::Row> BinRelWT_sdsl::get_column(Column column) const {
    assert(column < num_columns_);

    uint64_t num_relations_in_col = wt_.rank(wt_.size(), column);

    std::vector<Row> column_vec;
    for (uint64_t i = 0; i < num_relations_in_col; ++i) {
        uint64_t col_index_in_wt = wt_.select(i + 1, column);
        uint64_t col_index_in_del = delimiters_.select0(col_index_in_wt + 1);
        uint64_t row_num = delimiters_.rank1(col_index_in_del) - 1;
        column_vec.push_back(row_num);
    }
    return column_vec;
}

bool BinRelWT_sdsl::get(Row row, Column column) const {
    // Not counting i first 1s indicating the boundries.
    uint64_t first_string_index = delimiters_.select1(row + 1) - row;
    // Not counting i+1 first 1s indicating the boundries.
    uint64_t last_string_index = delimiters_.select1(row + 2) - row - 1;

    // Rank operation in sdsl::wt_in is exclusive.
    return first_string_index == 0
            ? wt_.rank(last_string_index, column)
            : wt_.rank(first_string_index, column)
                    != wt_.rank(last_string_index, column);
}

bool BinRelWT_sdsl::load(std::istream &in) {
    if (!in.good())
        return false;
    try {
        num_columns_ = load_number(in);
        wt_.load(in);
        return delimiters_.load(in);
    } catch (...) {
        return false;
    }
}

void BinRelWT_sdsl::serialize(std::ostream &out) const {
    if (!out.good())
        throw std::ofstream::failure("Bad stream");

    serialize_number(out, num_columns_);
    wt_.serialize(out);
    delimiters_.serialize(out);
}

uint64_t BinRelWT_sdsl::num_relations() const {
    return wt_.size();
}
