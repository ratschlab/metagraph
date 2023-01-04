#ifndef __INT_ROW_DISK_HPP__
#define __INT_ROW_DISK_HPP__

#include <vector>

#include "annotation/int_matrix/base/int_matrix.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"
#include "annotation/int_matrix/row_disk/disk_utils.h"

namespace mtg {
namespace annot {
namespace matrix {

//mkokot_TODO: consider merging it to a single class template with RowDisk class
// the idea is that the classes will be almost the same
// the only difference will be that there is some additional data stored
// for now for simplicity I just use different class
// and also different base class
class IntRowDisk : public IntMatrix {

public:
    IntRowDisk(size_t RA_ivbuffer_size = 16'384) {
        buffer_params_.buff_size = RA_ivbuffer_size/8;
        if (buffer_params_.buff_size < 8)
            buffer_params_.buff_size = 8;
    }

    uint64_t num_columns() const override { return num_columns_; }
    uint64_t num_rows() const override { return num_rows_; }

    bool get(Row row, Column column) const override {
        return View(boundary_, buffer_params_.filename,
                    buffer_params_.offset, buffer_params_.buff_size, bits_for_col_id_, bits_for_value_)
                .get(row, column);
    }

    SetBitPositions get_row(Row row) const override {
        return View(boundary_, buffer_params_.filename,
                    buffer_params_.offset, buffer_params_.buff_size, bits_for_col_id_, bits_for_value_)
                .get_row(row);
    }

    std::vector<Row> get_column(Column column) const override {
        mtg::common::logger->warn(
                "get_column is extremely inefficient operation, consider using "
                "column-major format");
        return View(boundary_, buffer_params_.filename, buffer_params_.offset,
                    buffer_params_.buff_size, bits_for_col_id_, bits_for_value_)
                .get_column(num_rows(), column);
    }

    std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const override {
        return View(boundary_, buffer_params_.filename,
                    buffer_params_.offset, buffer_params_.buff_size, bits_for_col_id_, bits_for_value_)
                .get_rows(rows);
    }

    RowValues get_row_values(Row row) const override {
        return View(boundary_, buffer_params_.filename,
                    buffer_params_.offset, buffer_params_.buff_size, bits_for_col_id_, bits_for_value_)
                .get_row_values(row);
    }

    std::vector<RowValues>
    get_row_values(const std::vector<Row> &rows) const override {
        return View(boundary_, buffer_params_.filename,
                    buffer_params_.offset, buffer_params_.buff_size, bits_for_col_id_, bits_for_value_)
                .get_row_values(rows);
    }

    bool load(std::istream &f) override;

    void serialize(std::ostream &/*out*/) const override {
        throw std::runtime_error("Not implemented");
        //mkokot_TODO: implement
    }

    // number of ones in the matrix
    uint64_t num_relations() const override {
        return num_set_bits_;
    }


    static void
    serialize(const std::function<void(std::function<void(const RowValues &)>)> &call_rows_with_values,
              const std::string &filename,
              uint64_t num_cols,
              uint64_t num_set_bits,
              uint64_t num_rows,
              uint8_t int_width_for_counts);

  private:

    // For the multithreading to work properly the idea is to open int_vector_buffer<> per
    // each public method call that needs to work with int_vector_buffer<> The actual
    // implementation is wrapped in the class below to assure it will not call any methods
    // of RowDisk class which could result in multiple opens of int_vector_buffer
    class View {
        std::ifstream open_and_set_pos(const std::string filename, size_t offset) {
            std::ifstream in(filename, std::ios::binary);
            if (!in)
                throw std::ofstream::failure("Cannot open file " + filename);
            in.seekg(offset, std::ios::beg);
            return in;
        }
      public:
        View(const bit_vector_small &boundary,
             const std::string &filename,
             uint64_t offset,
             uint64_t buff_size,
             uint64_t bits_for_col_id,
             uint64_t bits_for_value)
            : boundary_(boundary),
              in_(open_and_set_pos(filename, offset)),
              set_bits_(in_, buff_size),
              bits_for_col_id_(bits_for_col_id),
              bits_for_value_(bits_for_value) {}

        bool get(Row row, Column column) const;
        BinaryMatrix::SetBitPositions get_row(Row row) const;
        std::vector<SetBitPositions> get_rows(const std::vector<Row> &row_ids) const;
        std::vector<BinaryMatrix::Row> get_column(uint64_t num_rows, Column column) const;

        RowValues get_row_values(Row row) const;
        std::vector<RowValues> get_row_values(const std::vector<Row> &row_ids) const;

      private:
        const bit_vector_small &boundary_;
        //sdsl::int_vector_buffer<> set_bits_;
        std::ifstream in_;
        mutable DiskRandomReader set_bits_;
        uint64_t bits_for_col_id_;
        uint64_t bits_for_value_;
    };

    struct int_vector_buffer_params {
        std::string filename;
        uint64_t offset;
        uint64_t buff_size = 16'384/8; //in uint64_t
    };
    int_vector_buffer_params buffer_params_;

    bit_vector_small boundary_;
    uint64_t num_columns_ = 0;
    uint64_t num_set_bits_ = 0;
    uint64_t bits_for_col_id_ = 0;
    uint64_t bits_for_value_ = 0;
    uint64_t num_rows_ = 0;

    size_t iv_size_on_disk_ = 0; // for non-static serialization
};

} // namespace matrix
} // namespace annot
} // namespace mtg


#endif // __INT_ROW_DISK_HPP__