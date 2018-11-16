#include "rainbowfish.hpp"

#include <algorithm>
#include <functional>
#include <tsl/hopscotch_map.h>

#include "utils.hpp"
#include "serialization.hpp"


Rainbowfish::Rainbowfish(const std::function<void(RowCallback)> &call_rows,
                         uint64_t num_columns,
                         uint64_t buffer_size)
      : num_columns_(num_columns),
        buffer_size_(buffer_size) {
    assert(buffer_size_);

    using IndexVectorMap = tsl::hopscotch_map<SmallVector,
                                              uint64_t,
                                              SmallVectorHash>;
    IndexVectorMap vector_counter;
    IndexVectorMap vector_coder;
    uint64_t coded_rows_size = 0;
    uint64_t coded_rows_set_bits = 0;
    uint64_t counter_num_set_bits = 0;

    auto flush = [&](const auto &callback) {
        if (!vector_counter.size())
            return;

        std::vector<SmallVector> vectors;
        vectors.reserve(vector_counter.size());

        std::vector<std::pair<uint64_t, uint64_t>> index_counts;
        index_counts.reserve(vector_counter.size());

        for (auto it = vector_counter.begin(); it != vector_counter.end(); ++it) {
            index_counts.emplace_back(index_counts.size(), it->second);
            vectors.push_back(std::move(it->first));
        }
        vector_counter.clear();
        counter_num_set_bits = 0;

        // sort in the order of decreasing counts
        std::sort(index_counts.begin(), index_counts.end(),
            [](const auto &first, const auto &second) {
                return first.second > second.second;
            }
        );

        std::vector<Column> v;
        for (const auto &index_count : index_counts) {
            v.assign(vectors[index_count.first].begin(),
                     vectors[index_count.first].end());

            callback(v);

            v.clear();

            uint64_t code = vector_coder.size();
            coded_rows_size += utils::code_length(code) * index_count.second;
            coded_rows_set_bits += sdsl::bits::cnt(code) * index_count.second;

            assert(vector_coder.find(vectors[index_count.first]) == vector_coder.end());
            vector_coder[std::move(vectors[index_count.first])] = code;
        }
    };

    uint64_t rows = 0;
    SmallVector row_indices_small;
    call_rows([&](const utils::SetBitPositions &row_indices) {
        ++rows;
        row_indices_small.assign(row_indices.begin(), row_indices.end());
        if (vector_coder.size()) {
            auto code_find = vector_coder.find(row_indices_small);
            if (code_find != vector_coder.end()) {
                coded_rows_size += utils::code_length(code_find.value());
                coded_rows_set_bits += sdsl::bits::cnt(code_find.value());
                return;
            }
        }

        auto counter_find = vector_counter.try_emplace(row_indices_small, 1);
        if (counter_find.second) {
            counter_num_set_bits += row_indices_small.size();
        } else {
            counter_find.first.value()++;
        }

        if (vector_counter.size() == buffer_size_) {
            assert(counter_num_set_bits <= num_columns_ * buffer_size_);
            reduced_matrix_.emplace_back(
                new ReducedMatrixType(flush,
                                      num_columns_,
                                      buffer_size_,
                                      counter_num_set_bits)
            );
        }
    });

    uint64_t remaining_rows = vector_counter.size();
    assert(counter_num_set_bits <= num_columns_ * remaining_rows);
    reduced_matrix_.emplace_back(
        new ReducedMatrixType(flush,
                              num_columns_,
                              remaining_rows,
                              counter_num_set_bits)
    );

    sdsl::bit_vector row_builder(coded_rows_size);
    sdsl::bit_vector delimiter_vector(coded_rows_size, false);

    uint64_t pos = 0;
    row_indices_small.clear();
    call_rows([&](const auto &row_indices) {
        row_indices_small.assign(row_indices.begin(), row_indices.end());

        num_relations_ += row_indices_small.size();

        uint64_t code = vector_coder[row_indices_small];
        uint64_t code_length = utils::code_length(code);

        assert(code % buffer_size_
            < reduced_matrix_[code / buffer_size_]->num_rows());

        delimiter_vector[pos + code_length - 1] = true;

        row_builder.set_int(pos, code, code_length);
        pos += code_length;
    });
    vector_coder.clear();

    assert(pos == row_builder.size());

    row_code_delimiters_ = bit_vector_rrr<>(std::move(delimiter_vector));
    row_codes_ = decltype(row_codes_)(std::move(row_builder));
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
    uint64_t begin = row ? row_code_delimiters_.select1(row) + 1: 0;
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
std::vector<Rainbowfish::Column> Rainbowfish::get_row(Row row) const {
    uint64_t code = get_code(row);
    return reduced_matrix_[code / buffer_size_]->get_row(code % buffer_size_);
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
        row_codes_.load(in);
        row_code_delimiters_.load(in);

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
