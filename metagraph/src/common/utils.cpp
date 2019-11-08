#include "utils.hpp"

#include <cstdio>
#include <fstream>
#include <algorithm>
#include <filesystem>

#include "serialization.hpp"


namespace utils {

static bool VERBOSE = false;

bool get_verbose() { return VERBOSE; }
void set_verbose(bool verbose) { VERBOSE = verbose; }


RowsFromColumnsTransformer
::RowsFromColumnsTransformer(uint64_t num_rows,
                             const std::vector<std::string> &files)
      : num_rows_(num_rows) {
    streams_.reserve(files.size());
    for (const auto &file : files) {
        // initialize stream
        streams_.emplace_back(new VectorFileStream(file));
    }

    init_heap();
}

template <typename BitVectorPtr>
RowsFromColumnsTransformer
::RowsFromColumnsTransformer(const std::vector<BitVectorPtr> &columns) {
    if (!columns.size()) {
        num_rows_ = 0;
        return;
    }

    num_rows_ = columns.at(0)->size();
    for (const auto &column : columns) {
        if (column->size() != num_rows_)
            throw std::runtime_error("Error: columns are not of the same size");
    }

    streams_.reserve(columns.size());
    for (const auto &col_ptr : columns) {
        streams_.emplace_back(new VectorBitStream(*col_ptr));
    }

    init_heap();
}

template
RowsFromColumnsTransformer
::RowsFromColumnsTransformer(const std::vector<const bit_vector_small*> &);

template
RowsFromColumnsTransformer
::RowsFromColumnsTransformer(const std::vector<std::unique_ptr<bit_vector_sd>> &);

template
RowsFromColumnsTransformer
::RowsFromColumnsTransformer(const std::vector<std::unique_ptr<bit_vector_small>> &);

template
RowsFromColumnsTransformer
::RowsFromColumnsTransformer(const std::vector<std::unique_ptr<bit_vector>> &);

template
RowsFromColumnsTransformer
::RowsFromColumnsTransformer(const std::vector<std::shared_ptr<bit_vector_small>> &);


RowsFromColumnsTransformer
::RowsFromColumnsTransformer(const bit_vector_small &columns_concatenated,
                             uint64_t column_size) {
    assert(columns_concatenated.size() % column_size == 0);

    if (!columns_concatenated.size()) {
        num_rows_ = 0;
        return;
    }

    num_rows_ = column_size;
    uint64_t num_columns = columns_concatenated.size() / column_size;
    streams_.reserve(num_columns);
    for (uint64_t i = 0; i < num_columns; ++i) {
        streams_.emplace_back(new VectorBitStream(columns_concatenated,
                                                  column_size * i,
                                                  column_size * (i + 1)));
    }

    init_heap();
}

void RowsFromColumnsTransformer::init_heap() {
    for (uint64_t column = 0; column < streams_.size(); ++column) {
        assert(streams_[column].get());

        num_set_bits_left_ += streams_[column]->values_left();

        if (streams_[column]->values_left()) {
            index_heap_.emplace(kmer_label_pair {
                streams_[column]->next_value(),
                column
            });
        }
    }
}

void RowsFromColumnsTransformer::call_next(ValueCallback callback) {
    assert(values_left() && index_heap_.size());

    auto index = std::move(index_heap_.top());

    index_heap_.pop();
    num_set_bits_left_--;

    // TODO: Important!
    // Load multiple values from the stream to improve cache efficiency
    if (streams_[index.col_id]->values_left())
        index_heap_.emplace(kmer_label_pair {
            streams_[index.col_id]->next_value(),
            index.col_id
        });

    assert(index_heap_.size() || !values_left());
    callback(index.row_id, index.col_id);
}

std::tuple<uint64_t, uint64_t> RowsFromColumnsIterator::next_set_bit() {
    uint64_t row;
    uint64_t column;
    transformer_->call_next([&](uint64_t row_, uint64_t column_) {
        row = row_;
        column = column_;
    });
    return std::make_tuple(row, column);
}

std::vector<uint64_t> RowsFromColumnsIterator::next_row() {
    std::vector<uint64_t> indices;

    if (i_ > 0 && (row_ == i_)) {
        indices.push_back(column_);
    }

    if (!values_left() || row_ > i_) {
        i_++;
        return indices;
    }

    while (values_left()) {
        std::tie(row_, column_) = next_set_bit();
        if (row_ != i_)
            break;
        indices.push_back(column_);
    }
    i_++;
    return indices;
}

void call_rows(const std::function<void(const std::vector<uint64_t>&)> &callback,
               RowsFromColumnsTransformer&& transformer) {
    uint64_t cur_row = 0;
    std::vector<uint64_t> indices;

    while (transformer.values_left()) {
        transformer.call_next([&](uint64_t row, uint64_t column) {
            while (cur_row < row) {
                callback(indices);
                indices.clear();
                cur_row++;
            }
            indices.push_back(column);
        });
    }

    while (cur_row++ < transformer.rows()) {
        callback(indices);
        indices.clear();
    }

    assert(!transformer.values_left());
}

template <class BitVectorType>
std::vector<std::unique_ptr<bit_vector>>
transpose(const std::vector<std::unique_ptr<bit_vector>> &matrix) {
    std::vector<std::unique_ptr<bit_vector>> transposed;
    if (!matrix.size())
        return transposed;

    uint64_t num_rows = matrix.size();
    uint64_t num_columns = matrix[0]->size();
    transposed.reserve(num_columns);

    utils::call_rows(
        [&](const std::vector<uint64_t> &column_indices) {
            sdsl::bit_vector bv(num_rows, false);
            for (const auto &row_id : column_indices) {
                bv[row_id] = true;
            }
            transposed.emplace_back(new BitVectorType(std::move(bv)));
        },
        matrix
    );
    return transposed;
}

template
std::vector<std::unique_ptr<bit_vector>>
transpose<bit_vector_stat>(const std::vector<std::unique_ptr<bit_vector>> &matrix);

template
std::vector<std::unique_ptr<bit_vector>>
transpose<bit_vector_small>(const std::vector<std::unique_ptr<bit_vector>> &matrix);

template
std::vector<std::unique_ptr<bit_vector>>
transpose<bit_vector_dyn>(const std::vector<std::unique_ptr<bit_vector>> &matrix);


// indexes are distinct and sorted
sdsl::bit_vector subvector(const bit_vector &col,
                           const std::vector<uint64_t> &indexes) {
    assert(indexes.size() <= col.size());

    sdsl::bit_vector shrinked(indexes.size(), 0);

    uint64_t max_rank = col.num_set_bits();
    if (!max_rank)
        return shrinked;

    // the case of uncompressed vector
    if (dynamic_cast<const bit_vector_stat *>(&col)) {
        for (size_t j = 0; j < indexes.size(); ++j) {
            if (col[indexes[j]])
                shrinked[j] = true;
        }
        return shrinked;
    }

    uint64_t cur_rank = 1;
    uint64_t next_pos = col.select1(1);

    for (size_t j = 0; j < indexes.size(); ++j) {
        if (indexes[j] < next_pos)
            continue;

        if (indexes[j] == next_pos) {
            shrinked[j] = true;
            continue;
        }

        // indexes[j] > next_pos
        if (col[indexes[j]]) {
            shrinked[j] = true;
            continue;
        }

        // we found a zero, update next_pos
        cur_rank = col.rank1(indexes[j]) + 1;
        if (cur_rank > max_rank)
            return shrinked;

        next_pos = col.select1(cur_rank);
    }

    return shrinked;
}

std::vector<uint64_t> sample_indexes(uint64_t universe_size,
                                     uint64_t sample_size,
                                     std::mt19937 &gen) {
    if (!universe_size)
        return {};

    sample_size = std::min(universe_size, sample_size);

    std::vector<uint64_t> indexes;
    indexes.reserve(3 * sample_size);

    if (sample_size * 10 < universe_size) {
        std::uniform_int_distribution<uint64_t> dis(0, universe_size - 1);

        while (indexes.size() < sample_size) {
            indexes.clear();
            for (size_t i = 0; i < 1.5 * sample_size; ++i) {
                indexes.push_back(dis(gen));
            }
            std::sort(indexes.begin(), indexes.end());
            indexes.erase(std::unique(indexes.begin(), indexes.end()),
                          indexes.end());
        }
    } else {
        std::bernoulli_distribution dis(2.0 * sample_size / universe_size);

        while (indexes.size() < sample_size) {
            indexes.clear();
            for (size_t i = 0; i < universe_size; ++i) {
                if (dis(gen)) {
                    indexes.push_back(i);
                }
            }
        }
    }

    std::shuffle(indexes.begin(), indexes.end(), gen);

    return std::vector<uint64_t>(indexes.begin(), indexes.begin() + sample_size);
}

} // namespace utils


template <typename T>
std::set<T> convert_to_set(const std::vector<T> &vector) {
    return std::set<T>(vector.begin(), vector.end());
}

template
std::set<std::string>
convert_to_set(const std::vector<std::string> &vector);

template
std::set<uint64_t>
convert_to_set(const std::vector<uint64_t> &vector);


std::set<std::string> convert_to_set(const std::vector<std::string> &vector) {
    return convert_to_set<std::string>(vector);
}

std::set<std::pair<std::string, size_t>>
to_set(const std::vector<std::pair<std::string, size_t>> &vector) {
    return std::set<std::pair<std::string, size_t>>(vector.begin(), vector.end());
}
