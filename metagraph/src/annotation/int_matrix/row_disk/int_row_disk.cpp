#include "int_row_disk.hpp"

#include "common/logger.hpp"
#include "common/utils/file_utils.hpp"

namespace mtg {
namespace annot {
namespace matrix {

using mtg::common::logger;

bool IntRowDisk::View::get(Row row, Column column) const {
    SetBitPositions set_bits = get_row(row);
    SetBitPositions::iterator v = std::lower_bound(set_bits.begin(), set_bits.end(), column);
    return v != set_bits.end() && *v == column;
}

std::vector<IntMatrix::Row> IntRowDisk::View::get_column(uint64_t num_rows,
                                                         Column column) const {
    std::vector<Row> result;
    for (Row row = 0; row < num_rows; ++row) {
        if (get(row, column))
            result.push_back(row);
    }
    return result;
}

IntMatrix::SetBitPositions IntRowDisk::View::get_row(Row row) const {
    assert(boundary_[boundary_.size() - 1] == 1);
    uint64_t start_idx = row == 0 ? 0 : boundary_.select1(row) + 1;
    uint64_t end_idx = boundary_.next1(start_idx);

    SetBitPositions result;
    result.reserve(end_idx - start_idx);
    // In each row, the first value in `set_bits_` stores the first set bit,
    // and all next ones store deltas pos[i] - pos[i-1].
    uint64_t last = 0;    
    set_bits_.start_reading_at((start_idx - row) * (bits_for_col_id_ + bits_for_value_));
    Column col_id;
    Value val;
    for (uint64_t i = start_idx; i != end_idx; ++i) {        
        set_bits_.get(col_id, bits_for_col_id_);
        set_bits_.get(val, bits_for_value_); // skip
        result.push_back(last += col_id);
    }

    return result;
}

std::vector<IntMatrix::SetBitPositions>
IntRowDisk::View::get_rows(const std::vector<Row> &row_ids) const {
    std::vector<SetBitPositions> rows(row_ids.size());

    for (size_t i = 0; i < row_ids.size(); ++i) {
        rows[i] = get_row(row_ids[i]);
    }

    return rows;
}

IntMatrix::RowValues IntRowDisk::View::get_row_values(Row row) const {    
    assert(boundary_[boundary_.size() - 1] == 1);
    uint64_t start_idx = row == 0 ? 0 : boundary_.select1(row) + 1;
    uint64_t end_idx = boundary_.next1(start_idx);

    RowValues result;
    result.reserve(end_idx - start_idx);
    // In each row, the first value in `set_bits_` stores the first set bit,
    // and all next ones store deltas pos[i] - pos[i-1].
    uint64_t last = 0;
    set_bits_.start_reading_at((start_idx - row) * (bits_for_col_id_ + bits_for_value_));
    Column col_id;
    Value val;
    for (uint64_t i = start_idx; i != end_idx; ++i) {        
        set_bits_.get(col_id, bits_for_col_id_);
        set_bits_.get(val, bits_for_value_);         
        result.emplace_back(last += col_id, val);
    }

    return result;
}


std::vector<IntMatrix::RowValues>
IntRowDisk::View::get_row_values(const std::vector<Row> &row_ids) const {    
    std::vector<RowValues> rows_with_values(row_ids.size());

    for (size_t i = 0; i < row_ids.size(); ++i) {
        rows_with_values[i] = get_row_values(row_ids[i]);
    }
    
    return rows_with_values;
}

bool IntRowDisk::load(std::istream &f) {    
    auto _f = dynamic_cast<utils::NamedIfstream *>(&f);
    assert(_f);
    try {
        num_columns_ = load_number(f);
        auto boundary_start = load_number(f);

        num_set_bits_ = load_number(f);
        bits_for_col_id_ = load_number(f);
        bits_for_value_ = load_number(f);

        buffer_params_.filename = _f->get_name();
        buffer_params_.offset = f.tellg();

        assert(boundary_start >= buffer_params_.offset);
        iv_size_on_disk_ = boundary_start - buffer_params_.offset;

        f.seekg(boundary_start, std::ios_base::beg);

        boundary_.load(f);

        num_rows_ = boundary_.num_set_bits();

    } catch (...) {
        return false;
    }

    return true;
}

void IntRowDisk::serialize(
        const std::function<void(std::function<void(const RowValues &)>)> &call_rows_with_values,
        const std::string &filename,
        uint64_t num_cols,
        uint64_t num_set_bits,
        uint64_t num_rows,
        uint8_t int_width_for_counts) {
    // std::ios::ate needed because labels are serialized before
    // std::ios::in needed for tellp to return absolute file position
    std::ofstream outstream(filename, std::ios::binary | std::ios::ate | std::ios::in);

    if (!outstream.good())
        throw std::ofstream::failure("Cannot write to file " + filename);

    serialize_number(outstream, num_cols);

    auto boundary_start_pos = outstream.tellp();

    // boundary is built in memory
    // the boundary must be loaded in the IntRowDisk::load so its position must be known
    // binary format for this representation is
    // [boundary pos] (8B)
    // [serialized rows with values]
    // [boundary] (that are stored in the file at boundary pos)

    // write "empty" boundary start
    std::streampos boundary_start = 0;
    serialize_number(outstream, boundary_start);

    serialize_number(outstream, num_set_bits);

    uint8_t bits_for_col_id = sdsl::bits::hi(num_cols) + 1;
    uint8_t bits_for_value = int_width_for_counts;
    //uint8_t int_width = std::max(bits_for_col_id, bits_for_value);

    serialize_number(outstream, bits_for_col_id);
    serialize_number(outstream, bits_for_value);

    DiskWriter writer(outstream, 1024 * 1024);
    
    //  call_bits should be called for each argument that should be set (row boundaries)
    auto call_bits = [&](const std::function<void(uint64_t)> &call_bit) {
        uint64_t t = 0;        
        call_rows_with_values([&](const auto &row_with_values) {
            uint64_t last = 0;
            for (auto [col_id, val] : row_with_values) {                
                writer.add(col_id - last, bits_for_col_id);
                //outbuf.push_back(col_id - last);
                writer.add(val, bits_for_value);
                //outbuf.push_back(val);
                last = col_id;                
            }
            call_bit((t += row_with_values.size() + 1) - 1);
        });
        writer.flush();
    };
    
    bit_vector_small boundary(call_bits, num_rows + num_set_bits, num_rows);

    boundary_start = outstream.tellp();    
    boundary.serialize(outstream);
    outstream.seekp(boundary_start_pos, std::ios_base::beg);
    serialize_number(outstream, boundary_start);
}

} // namespace binmat
} // namespace annot
} // namespace mtg
