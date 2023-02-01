#ifndef __INT_COORD_DISK_HPP__
#define __INT_COORD_DISK_HPP__

#include <vector>

#include "annotation/int_matrix/base/int_matrix.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/disk_buffer.hpp"

namespace mtg {
namespace annot {
namespace matrix {

//mkokot_TODO: consider merging it to a single class template with RowDisk class
// the idea is that the classes will be almost the same
// the only difference will be that there is some additional data stored
// for now for simplicity I just use different class
// and also different base class
class CoordRowDisk : public MultiIntMatrix {
  public:
    CoordRowDisk(size_t RA_ivbuffer_size = 16'384) {
        buffer_params_.buff_size = std::max((size_t)8, RA_ivbuffer_size / 8);
    }

    uint64_t num_columns() const { return num_columns_; }
    uint64_t num_rows() const { return num_rows_; }

    bool get(Row i, Column j) const { return get_view().get(i, j); }
    SetBitPositions get_row(Row i) const { return get_view().get_row(i); }
    // FYI: `get_column` is very inefficient, consider using column-major formats
    std::vector<Row> get_column(Column j) const { return get_view().get_column(j); }
    std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const {
        return get_view().get_rows(rows);
    }

    RowValues get_row_values(Row i) const { return get_view().get_row_values(i); }
    std::vector<RowValues> get_row_values(const std::vector<Row> &rows) const {
        return get_view().get_row_values(rows);
    }

    // return total number of attributes in all tuples
    uint64_t num_attributes() const { return num_attributes_; }

    // return entries of the matrix -- where each entry is a set of integers
    RowTuples get_row_tuples(Row i) const { return get_view().get_row_tuples(i); }
    std::vector<RowTuples> get_row_tuples(const std::vector<Row> &rows) const {
        return get_view().get_row_tuples(rows);
    }

    bool load(std::istream &f);

    void serialize(std::ostream &/*out*/) const {
        throw std::runtime_error("Not implemented");
        //mkokot_TODO: implement
    }

    // number of ones in the matrix
    uint64_t num_relations() const { return num_set_bits_; }

    static void serialize(
            const std::string &filename,
            const std::function<void(std::function<void(const RowTuples &)>)> &call_rows,
            uint64_t num_cols,
            uint64_t num_rows,
            uint64_t num_set_bits,
            uint64_t num_values,
            uint64_t max_val,
            uint64_t max_tuple_size);

  private:
    // For the multithreading to work properly, we open int_vector_buffer<> in
    // a special View class that has an actual implementation of the method.
    class View {
      public:
        View(const bit_vector_small &boundary, //sd_vector or int_vector_buffer at the begining (instead of bit_vector_small)
             const std::string &filename,
             uint64_t offset,
             uint64_t buff_size,
             uint64_t bits_for_col_id,
             uint64_t bits_for_number_of_vals,
             uint64_t bits_for_single_value)
            : boundary_(boundary),
              in_(open_and_set_pos(filename, offset)),
              set_bits_(in_, buff_size),
              bits_for_col_id_(bits_for_col_id),
              bits_for_number_of_vals_(bits_for_number_of_vals),
              bits_for_single_value_(bits_for_single_value) {}

        bool get(Row i, Column j) const;
        BinaryMatrix::SetBitPositions get_row(Row i) const;
        std::vector<SetBitPositions> get_rows(const std::vector<Row> &row_ids) const;
        std::vector<BinaryMatrix::Row> get_column(Column j) const;

        RowValues get_row_values(Row i) const;
        std::vector<RowValues> get_row_values(const std::vector<Row> &row_ids) const;

        RowTuples get_row_tuples(Row i) const;

        std::vector<RowTuples> get_row_tuples(const std::vector<Row> &rows) const;

      private:
        std::ifstream open_and_set_pos(const std::string filename, size_t offset) {
            std::ifstream in(filename, std::ios::binary);
            if (!in)
                throw std::ofstream::failure("Cannot open file " + filename);
            in.seekg(offset, std::ios::beg);
            return in;
        }

        const bit_vector_small &boundary_;
        std::ifstream in_;
        // layout: [col_id|value]
        // `boundary_` puts 0 for each col-value pair and delimits rows with 1
        mutable common::DiskRandomReader set_bits_;
        uint64_t bits_for_col_id_;
        uint64_t bits_for_number_of_vals_;
        uint64_t bits_for_single_value_;
    };

    View get_view() const {
        return View(boundary_, buffer_params_.filename,
                    buffer_params_.offset, buffer_params_.buff_size, bits_for_col_id_,
                    bits_for_number_of_vals_, bits_for_single_value_);
    }

    struct {
        std::string filename;
        uint64_t offset;
        uint64_t buff_size = 16'384 / 8; // in uint64_t
    } buffer_params_;

    bit_vector_small boundary_;
    uint64_t num_columns_ = 0;
    uint64_t num_set_bits_ = 0;
    uint64_t num_attributes_ = 0;
    uint64_t bits_for_col_id_ = 0;
    uint64_t bits_for_number_of_vals_ = 0;
    uint64_t bits_for_single_value_ = 0;
    uint64_t num_rows_ = 0;
};

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __INT_COORD_DISK_HPP__
