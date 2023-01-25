#ifndef __ROW_DISK_HPP__
#define __ROW_DISK_HPP__

#include <string>
#include <vector>

#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"

namespace mtg {
namespace annot {
namespace binmat {

class RowDisk : public BinaryMatrix {
  public:
    RowDisk(size_t RA_ivbuffer_size = 16'384) {
        buffer_params_.buff_size = RA_ivbuffer_size;
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

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    // number of ones in the matrix
    uint64_t num_relations() const { return num_relations_; }

    static void serialize(const std::string &filename,
                          const std::function<void(binmat::BinaryMatrix::RowCallback)> &call_rows,
                          uint64_t num_cols,
                          uint64_t num_set_bits,
                          uint64_t num_rows);

    const bit_vector_small& get_boundary() const { return boundary_; }

  private:
    // For the multithreading to work properly, we open int_vector_buffer<> in
    // a special View class that has an actual implementation of the method.
    class View {
      public:
        View(const bit_vector_small &boundary,
             const std::string &filename,
             uint64_t offset,
             uint64_t buff_size)
            : boundary_(boundary),
              set_bits_(filename, std::ios::in, buff_size, 0, false, offset) {}

        bool get(Row i, Column j) const;
        BinaryMatrix::SetBitPositions get_row(Row i) const;
        std::vector<SetBitPositions> get_rows(const std::vector<Row> &row_ids) const;
        std::vector<BinaryMatrix::Row> get_column(Column j) const;

      private:
        const bit_vector_small &boundary_;
        sdsl::int_vector_buffer<> set_bits_;
    };

    View get_view() const {
        return View(boundary_, buffer_params_.filename,
                    buffer_params_.offset, buffer_params_.buff_size);
    }

    struct {
        std::string filename;
        uint64_t offset;
        uint64_t buff_size = 16'384;
    } buffer_params_;

    bit_vector_small boundary_;
    uint64_t num_columns_ = 0;
    uint64_t num_rows_ = 0;
    uint64_t num_relations_ = 0;

    size_t iv_size_on_disk_ = 0; // for non-static serialization
};

} // namespace binmat
} // namespace annot
} // namespace mtg

#endif // __ROW_DISK_HPP__
