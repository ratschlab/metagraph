#include "binary_matrix.hpp"

#include <ips4o.hpp>

#include "common/algorithms.hpp"
#include "common/vectors/bitmap.hpp"
#include "common/serialization.hpp"
#include "common/utils/template_utils.hpp"
#include "common/vector_set.hpp"
#include "common/hashers/hash.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"


namespace mtg {
namespace annot {
namespace matrix {

const size_t kRowBatchSize = 10'000;
using Column = RainbowMatrix::Column;


std::vector<BinaryMatrix::SetBitPositions>
BinaryMatrix::get_rows_dict(std::vector<Row> *rows, size_t num_threads) const {
    VectorSet<SetBitPositions, utils::VectorHash> unique_rows;

    std::vector<std::pair<Row, size_t>> row_to_index(rows->size());
    for (size_t i = 0; i < rows->size(); ++i) {
        row_to_index[i] = std::make_pair((*rows)[i], i);
    }

    // don't break the topological order for row-diff annotation
    if (!dynamic_cast<const IRowDiff *>(this)) {
        ips4o::parallel::sort(row_to_index.begin(), row_to_index.end(),
                              utils::LessFirst(), num_threads);
    }

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (uint64_t begin = 0; begin < row_to_index.size(); begin += kRowBatchSize) {
        const uint64_t end = std::min(begin + kRowBatchSize,
                                      static_cast<uint64_t>(row_to_index.size()));

        std::vector<Row> ids(end - begin);
        for (uint64_t i = begin; i < end; ++i) {
            ids[i - begin] = row_to_index[i].first;
        }

        auto batch = get_rows(ids);

        #pragma omp critical
        {
            for (uint64_t i = begin; i < end; ++i) {
                auto it = unique_rows.emplace(std::move(batch[i - begin])).first;
                (*rows)[row_to_index[i].second] = it - unique_rows.begin();
            }
        }
    }

    return to_vector(std::move(unique_rows));
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
    min_count = std::max(min_count, (size_t)1);

    auto ic = index_counts;

    // don't break the topological order for row-diff annotation
    if (!dynamic_cast<const IRowDiff *>(this))
        std::sort(ic.begin(), ic.end(), utils::LessFirst());

    auto call_bits = [&](const auto &callback) {
        // limit the RAM usage by avoiding querying all the rows at once
        for (uint64_t begin = 0; begin < ic.size(); begin += kRowBatchSize) {
            const uint64_t end = std::min(begin + kRowBatchSize,
                                          static_cast<uint64_t>(ic.size()));

            std::vector<Row> ids(end - begin);
            for (uint64_t i = begin; i < end; ++i) {
                ids[i - begin] = ic[i].first;
            }

            const auto &rows = get_rows(ids);

            for (uint64_t i = begin; i < end; ++i) {
                for (size_t j : rows[i - begin]) {
                    callback(j, ic[i].second);
                }
            }
        }
    };

    return utils::accumulate_counts(call_bits, num_columns(), min_count);
}


std::vector<RainbowMatrix::SetBitPositions>
RainbowMatrix::get_rows(const std::vector<Row> &rows) const {
    std::vector<Row> pointers = rows;
    auto distinct_rows = get_rows_dict(&pointers);

    std::vector<SetBitPositions> result(rows.size());
    for (size_t i = 0; i < pointers.size(); ++i) {
        result[i] = distinct_rows[pointers[i]];
    }

    return result;
}

std::vector<RainbowMatrix::SetBitPositions>
RainbowMatrix::get_rows_dict(std::vector<Row> *rows, size_t num_threads) const {
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

    if (num_threads <= 1)
        return codes_to_rows(codes);

    std::vector<SetBitPositions> unique_rows(codes.size());

    size_t batch_size = std::min(kRowBatchSize,
                                 (codes.size() + num_threads - 1) / num_threads);

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (size_t i = 0; i < codes.size(); i += batch_size) {
        std::vector<uint64_t> ids(codes.begin() + i,
                                  codes.begin() + std::min(i + batch_size, codes.size()));
        auto rows = codes_to_rows(ids);
        for (size_t j = 0; j < rows.size(); ++j) {
            unique_rows[i + j] = std::move(rows[j]);
        }
    }

    return unique_rows;
}

std::vector<std::pair<RainbowMatrix::Column, size_t /* count */>>
RainbowMatrix::sum_rows(const std::vector<std::pair<Row, size_t>> &index_counts,
                        size_t min_count) const {
    min_count = std::max(min_count, (size_t)1);

    std::vector<Row> row_ids;
    row_ids.reserve(index_counts.size());

    size_t total_sum_count = 0;
    for (auto [i, count] : index_counts) {
        row_ids.push_back(i);
        total_sum_count += count;
    }

    if (total_sum_count < min_count)
        return {};

    auto distinct_rows = get_rows_dict(&row_ids);

    Vector<size_t> weights(distinct_rows.size(), 0);
    for (size_t i = 0; i < index_counts.size(); ++i) {
        weights[row_ids[i]] += index_counts[i].second;
    }

    auto call_bits = [&](const auto &callback) {
        for (size_t i = 0; i < distinct_rows.size(); ++i) {
            for (Column j : distinct_rows[i]) {
                callback(j, weights[i]);
            }
        }
    };
    return utils::accumulate_counts(call_bits, num_columns(), min_count);
}


std::vector<BinaryMatrix::SetBitPositions>
RowMajor::get_rows(const std::vector<Row> &row_ids) const {
    std::vector<SetBitPositions> rows(row_ids.size());
    for (size_t i = 0; i < row_ids.size(); ++i) {
        rows[i] = get_row(row_ids[i]);
    }
    return rows;
}

} // namespace matrix
} // namespace annot
} // namespace mtg
