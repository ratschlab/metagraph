#include "coord_row_disk.hpp"

#include "common/logger.hpp"
#include "common/utils/file_utils.hpp"
#include "common/ifstream_with_name.hpp"


namespace mtg {
namespace annot {
namespace matrix {

bool CoordRowDisk::View::get(Row row, Column column) const {
    SetBitPositions set_bits = get_row(row);
    SetBitPositions::iterator v = std::lower_bound(set_bits.begin(), set_bits.end(), column);
    return v != set_bits.end() && *v == column;
}

std::vector<IntMatrix::Row> CoordRowDisk::View::get_column(uint64_t num_rows,
                                                         Column column) const {
    std::vector<Row> result;
    for (Row row = 0; row < num_rows; ++row) {
        if (get(row, column))
            result.push_back(row);
    }
    return result;
}

IntMatrix::SetBitPositions CoordRowDisk::View::get_row(Row row) const {
    assert(boundary_[boundary_.size() - 1] == 1);
    uint64_t start_idx = row == 0 ? 0 : boundary_.select1(row) + 1;
    uint64_t end_idx = boundary_.next1(start_idx);

    SetBitPositions result;
    
    Column col_id;
    
    size_t pos = start_idx;
    while (pos < end_idx) {
        set_bits_.start_reading_at(pos - row);
        set_bits_.get(col_id, bits_for_col_id_);
        pos += bits_for_col_id_;
        
        size_t num_values;
        set_bits_.get(num_values, bits_for_number_of_vals_);

        assert(num_values);

        pos += bits_for_number_of_vals_;
        pos += bits_for_single_value_ * num_values;
        result.emplace_back(col_id);
    }

    return result;
}

std::vector<IntMatrix::SetBitPositions>
CoordRowDisk::View::get_rows(const std::vector<Row> &row_ids) const {
    std::vector<SetBitPositions> rows(row_ids.size());

    for (size_t i = 0; i < row_ids.size(); ++i) {
        rows[i] = get_row(row_ids[i]);
    }

    return rows;
}

IntMatrix::RowValues CoordRowDisk::View::get_row_values(Row row) const {
    assert(boundary_[boundary_.size() - 1] == 1);
    uint64_t start_idx = row == 0 ? 0 : boundary_.select1(row) + 1;
    uint64_t end_idx = boundary_.next1(start_idx);

    RowValues result;    
    Column col_id;    
    size_t pos = start_idx;
    while (pos < end_idx) {

        set_bits_.start_reading_at(pos - row);

        set_bits_.get(col_id, bits_for_col_id_);
        
        pos += bits_for_col_id_;
        
        size_t num_values;
        set_bits_.get(num_values, bits_for_number_of_vals_);

        assert(num_values);

        pos += bits_for_number_of_vals_;
        pos += num_values * bits_for_single_value_;

        result.emplace_back(col_id, num_values);
    }

    return result;
}


std::vector<IntMatrix::RowValues>
CoordRowDisk::View::get_row_values(const std::vector<Row> &row_ids) const {
   std::vector<RowValues> rows_with_values(row_ids.size());

    for (size_t i = 0; i < row_ids.size(); ++i) {
        rows_with_values[i] = get_row_values(row_ids[i]);
    }
    
    return rows_with_values;
}

MultiIntMatrix::RowTuples CoordRowDisk::View::get_row_tuples(Row row) const {
    assert(boundary_[boundary_.size() - 1] == 1);
    uint64_t start_idx = row == 0 ? 0 : boundary_.select1(row) + 1;
    uint64_t end_idx = boundary_.next1(start_idx);

    RowTuples result;
    
    set_bits_.start_reading_at(start_idx - row);
    Column col_id;
    
    size_t pos = start_idx;
    while (pos < end_idx) {
        set_bits_.get(col_id, bits_for_col_id_);
        pos += bits_for_col_id_;
        
        size_t num_values;
        set_bits_.get(num_values, bits_for_number_of_vals_);

        assert(num_values);

        pos += bits_for_number_of_vals_;
        Tuple tuple;
        tuple.reserve(num_values);
        for(size_t i = 0 ; i < num_values ; ++i) {
            Value val;
            set_bits_.get(val, bits_for_single_value_);
            tuple.push_back(val);
            pos+=bits_for_single_value_;
        }

        result.emplace_back(col_id, std::move(tuple));
    }

    return result;
}
        
std::vector<MultiIntMatrix::RowTuples> CoordRowDisk::View::get_row_tuples(const std::vector<Row> &rows) const {
    std::vector<RowTuples> rows_with_tuples(rows.size());

    for (size_t i = 0; i < rows.size(); ++i) {
        rows_with_tuples[i] = get_row_tuples(rows[i]);
    }
    
    return rows_with_tuples;
}


bool CoordRowDisk::load(std::istream &f) {
    auto _f = dynamic_cast<mtg::common::IfstreamWithName *>(&f);
    assert(_f);
    try {
        num_columns_ = load_number(f);
        auto boundary_start = load_number(f);

        num_set_bits_ = load_number(f);
        bits_for_col_id_ = load_number(f);
        bits_for_number_of_vals_ = load_number(f);
        bits_for_single_value_ = load_number(f);

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

void CoordRowDisk::serialize(
        const std::function<void(std::function<void(const RowTuples &)>)> &write_row_with_tuples,
        const std::string & filename,
        uint64_t num_cols,
        uint64_t num_set_bits,
        uint64_t num_rows,
        uint64_t num_values) {

    // std::ios::ate needed because labels are serialized before
    // std::ios::in needed for tellp to return absolute file position
    std::ofstream outstream(filename, std::ios::binary | std::ios::ate | std::ios::in);

    if (!outstream.good())
        throw std::ofstream::failure("Cannot write to file " + filename);

    serialize_number(outstream, num_cols);

    auto boundary_start_pos = outstream.tellp();

    std::streampos boundary_start = 0;
    serialize_number(outstream, boundary_start);

    serialize_number(outstream, num_set_bits);

    uint8_t bits_for_col_id = sdsl::bits::hi(num_cols) + 1;
    uint8_t bits_for_number_of_vals = 32; //mkokot_TODO: maybe would be possible to calculate better, or store for each row etc.
    uint8_t bits_for_single_value = 32; //mkokot_TODO: maybe would be possible to calculate better, or store for each row etc.
    
    serialize_number(outstream, bits_for_col_id);
    serialize_number(outstream, bits_for_number_of_vals);
    serialize_number(outstream, bits_for_single_value);

    DiskWriter writer(outstream, 1024 * 1024);

    auto call_bits = [&](const std::function<void(uint64_t)> &call_bit) {
        uint64_t t = 0;        
        write_row_with_tuples([&](const RowTuples &row_tuples) {            
            size_t bits_used = 0;
            
            for (const auto &tuple : row_tuples) {                
                auto col_id = tuple.first;
                const auto& values = tuple.second;

                writer.add(col_id, bits_for_col_id);
                bits_used += bits_for_col_id;

                writer.add(values.size(), bits_for_number_of_vals);
                bits_used += bits_for_number_of_vals;

                for (auto v : values) {
                    writer.add(v, bits_for_single_value);
                    bits_used += bits_for_single_value;
                }                
            }
            auto boundary_val = (t += bits_used + 1) - 1;            
            call_bit(boundary_val);            
        });
        writer.flush();
    };

    uint64_t boundary_capacity = num_set_bits * (bits_for_col_id + bits_for_number_of_vals) + num_values * bits_for_single_value + num_rows;

    bit_vector_small boundary(call_bits, boundary_capacity, num_rows);

    boundary_start = outstream.tellp();
    boundary.serialize(outstream);
    outstream.seekp(boundary_start_pos, std::ios_base::beg);
    serialize_number(outstream, boundary_start);    
}

} // namespace binmat
} // namespace annot
} // namespace mtg
