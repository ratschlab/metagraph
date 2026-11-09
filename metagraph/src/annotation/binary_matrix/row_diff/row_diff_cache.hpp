#pragma once

#include <cstdint>
#include <utility>

#include <tsl/hopscotch_map.h>

#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"


namespace mtg {
namespace annot {
namespace matrix {

/**
 * Decoded row-diff rows, kept across the reads of one reader. The row-diff paths of
 * neighbouring graph nodes share most of their rows, so a reader that walks the graph,
 * reading a few rows per step, decodes each path once instead of at every step.
 *
 * Bounded by an estimate of the bytes held: rows go into a new generation, and when it
 * holds half the bound it replaces the old one, which is dropped.
 */
template <class RowType>
class RowDiffCache {
  public:
    using Row = BinaryMatrix::Row;

    explicit RowDiffCache(uint64_t max_bytes = 64 << 20) : max_bytes_(max_bytes) {}

    const RowType* find(Row row) const {
        for (const Map &gen : gen_) {
            if (auto it = gen.find(row); it != gen.end())
                return &it->second;
        }
        return nullptr;
    }

    void insert(Row row, RowType&& value) {
        if (find(row))
            return;
        bytes_ += bytes_of(value);
        gen_[0].emplace(row, std::move(value));
        if (bytes_ > max_bytes_ / 2) {
            gen_[1].swap(gen_[0]);
            Map().swap(gen_[0]);
            bytes_ = 0;
        }
    }

    // queried rows whose path stopped at a kept row
    uint64_t hits = 0;

  private:
    using Map = tsl::hopscotch_map<Row, RowType>;

    static uint64_t bytes_of(const BinaryMatrix::SetBitPositions &row) {
        return 64 + row.size() * sizeof(row[0]);
    }
    static uint64_t bytes_of(const MultiIntMatrix::RowTuples &row) {
        uint64_t bytes = 64;
        for (const auto &[column, tuple] : row) {
            bytes += sizeof(row[0]) + tuple.size() * sizeof(tuple[0]);
        }
        return bytes;
    }

    Map gen_[2];
    uint64_t bytes_ = 0;   // of gen_[0]
    uint64_t max_bytes_;
};

} // namespace matrix
} // namespace annot
} // namespace mtg
