#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"

namespace mtg {
namespace annot {
namespace binmat {


class RowSparseDisk : public BinaryMatrix {
  public:
    RowSparseDisk(size_t RA_ivbuffer_size = 16'384) {
        buffer_params_.buff_size = RA_ivbuffer_size;
    }

    uint64_t num_columns() const override { return num_columns_; }
    uint64_t num_rows() const override { return num_rows_; }

    bool get(Row row, Column column) const override {
        return View(boundary_, buffer_params_.filename,
                    buffer_params_.offset, buffer_params_.buff_size)
                .get(row, column);
    }
    SetBitPositions get_row(Row row) const override {
        return View(boundary_, buffer_params_.filename,
                    buffer_params_.offset, buffer_params_.buff_size)
                .get_row(row);
    }
    std::vector<Row> get_column(Column column) const override {
        mtg::common::logger->warn(
                "get_column is extremely inefficient operation, consider using "
                "column-major format");
        return View(boundary_, buffer_params_.filename, buffer_params_.offset,
                    buffer_params_.buff_size)
                .get_column(num_rows(), column);
    }

    std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const override {
        return View(boundary_, buffer_params_.filename,
                    buffer_params_.offset, buffer_params_.buff_size)
                .get_rows(rows);
    }

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    // number of ones in the matrix
    uint64_t num_relations() const override {
        return View(boundary_, buffer_params_.filename,
                    buffer_params_.offset, buffer_params_.buff_size)
                .num_relations_impl();
    }

    static void serialize(const std::function<void(binmat::BinaryMatrix::RowCallback)> &call_rows,
                          const std::string &filename,
                          uint64_t num_cols,
                          uint64_t num_set_bits,
                          uint64_t num_rows);

    void set_buff_size(uint64_t buff_size) {
        buffer_params_.buff_size = buff_size;
    }

    const bit_vector_small get_boundary() const { return boundary_; }

  private:
    // For the multithreading to work properly the idea is to open int_vector_buffer<> per
    // each public method call that needs to work with int_vector_buffer<> The actual
    // implementation is wrapped in the class below to assure it will not call any methods
    // of RowSparseDisk class which could result in multiple opens of int_vector_buffer
    class View {
      public:
        View(const bit_vector_small &boundary,
             const std::string &filename,
             uint64_t offset,
             uint64_t buff_size)
            : boundary_(boundary),
              set_bits_(filename, std::ios::in, buff_size, 0, false, offset) {}

        uint64_t num_relations_impl() const { return set_bits_.size(); }
        bool get(Row row, Column column) const;
        BinaryMatrix::SetBitPositions get_row(Row row) const;
        std::vector<SetBitPositions> get_rows(const std::vector<Row> &row_ids) const;
        std::vector<BinaryMatrix::Row> get_column(uint64_t num_rows, Column column) const;

      private:
        const bit_vector_small &boundary_;
        sdsl::int_vector_buffer<> set_bits_;
    };

    struct int_vector_buffer_params {
        std::string filename;
        uint64_t offset;
        uint64_t buff_size = 10'000;
    };
    int_vector_buffer_params buffer_params_;

    bit_vector_small boundary_;
    uint64_t num_columns_ = 0;
    uint64_t num_rows_ = 0;

    size_t iv_size_on_disk_ = 0; // for non-static serialization
};


} // namespace binmat
} // namespace annot
} // namespace mtg
